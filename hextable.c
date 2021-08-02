/*  File: hextable.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) 1993-2021 Richard Durbin, Wellcome Sanger Institute
 *  License: available under the MIT license as in the accompanying LICENSE file 
 *-------------------------------------------------------------------
 * Description: makes hexamer tables for use by hexamer
                uses stats relative to composition only
 * Exported functions: main()
 * HISTORY:
 * Last edited: Aug  2 23:50 2021 (rd109)
 * * Aug  2 22:59 2021 (rd109): removed all acedb code in this standalone version
 * Created: Sun Aug 27 16:08:28 1995 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "readseq.h"
#include <stdarg.h>
#include <math.h>
#include <getopt.h>

int   nHex, hex[4096] ;
int   nHex2, hex2[4096] ;
int   codon[64] ;
float tab[4096] ;

int nseq = 0 ;
char *seqs[1000] ;
char *ids[1000] ;
int lens[1000] ;

int isCoding = 1 ;		/* default coding for hexExon */

void die (char *format, ...)
{
  va_list args ;

  if (strcmp (format, "usage"))
    { va_start (args, format) ;
      fprintf (stderr, "FATAL ERROR: ") ;
      vfprintf (stderr, format, args) ;
      fprintf (stderr, "\n") ;
      va_end (args) ;
    }
  fprintf (stderr, "Usage: hextable [-o ofile] [-2 file2] [-s sfile] file1\n") ;
  fprintf (stderr, "  all files are DNA fasta files\n") ;
  fprintf (stderr, "  -o <file>  output file\n") ;
  fprintf (stderr, "  -2 <file2> calculate stats by LLratio to file2\n") ;
  fprintf (stderr, "  -s <sfile> evaluates stats on sfile, not file1\n") ;
  fprintf (stdout, "  -n         flag for noncoding (no triplet frame)\n") ;

  exit (-1) ;
}

void information (int size, int *a)
{
  int i, max = 1 << (2*size) ;
  float x, I = 0 ;

  for (i = 0 ; i < max ; ++i) if (a[i])
    { x = a[i] / (float)nHex ;
      I -= x * log(x) ;
    }
  printf ("%.3f bits per base in %dmers\n", 
	  I / (size*log(2.0)), size) ;
}

void hexTableComposition (void)
{
  int i, j, k, index, nbad = 0, nstop = 0 ;
  float x, S = 0, E = 0, min = 0, max = 0 ;
  int comp[4096], compN[10000] ;
  float compSum[10000] ;
  float weight, order1 ;

			/* hexamers conditional on composition */

  memset (compSum, 0, 10000*sizeof(float)) ;
  memset (compN, 0, 10000*sizeof(int)) ;
  for (i = 0 ; i < 4096 ; ++i)
    { index = i ; k = 0 ;
      for (j = 6 ; j-- ; index >>= 2)
	switch (index & 0x3)
	  { 
	  case 0: k += 1 ;
	  case 1: k += 10 ;
	  case 2: k += 100 ;
	  case 3: k += 1000 ;
	  }
      comp[i] = k ;
      order1 = 1.0 ;		/* prob of hexamer i from 1st order model */
      hex[i] /= order1 ;	/* correct for 1st order biases */
      compSum[k] += hex[i] ;
      ++compN[k] ;
    }

  for (k = 0 ; k < 10000 ; ++k)
    if (compN[k])
      compSum[k] /= compN[k] ;

  for (i = 0 ; i < 4096 ; ++i)
    { j = i & 0x3f ; 
      if (isCoding && (j == 48 || j == 50 || j == 56))
	{ tab[i] = -100.0 ; ++nstop ; continue ; }
      j = i >> 6 ; 
      if (isCoding && (j == 48 || j == 50 || j == 56))
	{ tab[i] = -100.0 ; ++nstop ; continue ; }
      if (hex[i])
	{ k = comp[i] ;
	  x = hex[i] / compSum[k] ;
	  tab[i] = log(x) / log(2.0) ;
	  S += tab[i] * compSum[k] / nHex ;
	  E += tab[i] * (float)hex[i] / nHex ;
	  if (tab[i] > max) max = tab[i] ;
	  if (tab[i] < min) min = tab[i] ;
	}
      else
	{ tab[i] = -5.0 ;
	  ++nbad ;
	}
    }

  printf ("Hex table: %6.3f bits per triplet in coding\n"
	  "           %6.3f bits per triplet in scrambled coding\n",
	  0.5*E, 0.5*S) ;
  printf ("           min = %.2f, max = %.2f\n", min, max) ;
  if (isCoding)
    printf ("           %d stops, %d missing\n", nstop, nbad) ;
  else
    printf ("           %d missing\n", nbad) ;
}

void hexLikelihoodRatio ()
{ 
  int i, j ;
  float x, S = 0, E = 0, min = 0, max = 0 ;
  float rat = nHex2 / (float) nHex ;

  for (i = 0 ; i < 4096 ; ++i)
    { j = i & 0x3f ; 
      if (isCoding && (j == 48 || j == 50 || j == 56))
	{ tab[i] = -100.0 ; continue ; }
      j = i >> 6 ; 
      if (isCoding && (j == 48 || j == 50 || j == 56))
	{ tab[i] = -100.0 ; continue ; }
      x = rat * hex[i] / hex2[i] ;
      tab[i] = log (x) / log (2.0) ;
      S += tab[i] * hex[i] / (float) nHex ;
      E += tab[i] * hex2[i] / (float) nHex2 ;
      if (tab[i] > max) max = tab[i] ;
      if (tab[i] < min) min = tab[i] ;
    }
  printf ("LLR table: %6.3f bits per triplet in coding\n"
	  "           %6.3f bits per triplet in non-coding\n",
	  0.5*S, 0.5*E) ;
  printf ("           min = %.2f, max = %.2f\n", min, max) ;
}

float scoreSeq (int iseq, int debug)
{
  char *s ;
  int i, len, index ;
  float score = 0 ;

  s = seqs[iseq] ; len = lens[iseq] ;
  if (debug)
    { for (i = 0 ; i < len ; i += 3)
	printf ("  %c%c%c", index2char[s[i]], index2char[s[i+1]], index2char[s[i+2]]) ;
      printf ("\n     ") ;
    }
  index = (s[0] << 4) + (s[1] << 2) + s[2] ; s += 3 ;
  for (i = 3 ; i < len-3 ; i += 3, s += 3)
    { index &= 0x3f ;
      index <<= 6 ;
      index += (s[0] << 4) + (s[1] << 2) + s[2] ;
      score += tab[index] ;
      if (debug)
	{ if (tab[index] > -10)
	    printf (" %4.1f", tab[index]) ;
	  else
	    printf (" -100") ;
	}
    }

  if (debug)
    printf (" = %.2f\n", score) ;
      
  return score ;
}

void saveTable (char *name)
{
  FILE *fil ;
  int i, j, n = 0 ;

  if (!(fil = fopen (name, "w")))
    die ("Can't open output file %s", name) ;

  for ( i = 0 ; i < 256 ; ++i)
    { for (j = 0 ; j < 16 ; ++j, ++n)
	if (tab[n] > -10.0)
	  fprintf (fil, " %7.3f", tab[n]) ;
	else
	  fprintf (fil, " %7.1f", tab[n]) ;
      fprintf (fil, "\n") ;
    }

  fclose (fil) ;
}

void scoreSeqs (void)
{ 
  int iseq, i, nNeg = 0 ;
  float score, min = 0, max = 0, sumScore = 0 ;
  static int sHist[200] ;

  for (iseq = 0 ; iseq < nseq ; ++iseq)
    { score = scoreSeq (iseq, 0) ;
      sumScore += score ;
      if (score > max) max = score ;
      if (score < min) min = score ;
      if (score < 0) 
	{ nNeg++ ;
/*
	  printf ("%20s  %6.2f\n", ids[iseq], score) ;
	  if (score < -5)
	    { scoreSeq (iseq, 1) ;
	      ++seqs[iseq] ; --lens[iseq] ; scoreSeq (iseq, 1) ;
	      ++seqs[iseq] ; --lens[iseq] ; scoreSeq (iseq, 1) ;
	    }
*/
	}
      if (score < -1000) score = -1000 ;
      if (score > 999) score = 999 ;
      ++sHist[(int)(score+1000) / 10] ;
    }
  printf ("%d scores - average %.2f, max %.2f, min %.2f\n",
	  nseq, sumScore/nHex, max, min) ;
  printf ("            - %d less than 0\n", nNeg) ;
  for (i = 0 ; i < 200 ; ++i) if (sHist[i])
    printf ("  %3d :  %d\n", (i-100)*10, sHist[i]) ;
}

int main (int argc, char **argv)
{ 
  FILE *fil ;
  char *seq, *id, *s ;
  char *file1, *ofile = 0, *sfile = 0, *file2 = 0 ;
  int i, n, len, index ;

  while ((n = getopt (argc, argv, "o:s:2:n")) != -1)
    switch (n)
      {
      case 'o': ofile = optarg ; break ;
      case 's': sfile = optarg ; break ;
      case '2': file2 = optarg ; break ;
      case 'n': isCoding = 0 ; break ;
      default: die ("usage") ;
      }
  if (argc - optind != 1)
    die ("usage") ;
  file1 = argv[optind] ;

  dna2indexConv['n'] = dna2indexConv['N'] = -2 ;

  if (!(fil = fopen (file1, "r")))
    die ("Failed to open fasta file %s", file1) ;

				/* Dirichlet prior */
  for (i = 0 ; i < 4096 ; ++i)
    hex[i] = 1 ;
  nHex = 4096 ;

  while (readSequence (fil, dna2indexConv, &seq, &id, 0, &len))
    { s = seq ;
      seqs[nseq] = seq ; ids[nseq] = id ; lens[nseq] = len ; 
      if (++nseq > 1000)
	die ("More than 1000 sequences - edit and recompile") ;
      index = (s[0] << 4) + (s[1] << 2) + s[2] ; s += 3 ;
      for (i = 3 ; i < len-3 ; i += 3, s += 3)
	{ index &= 0x3f ;
	  ++codon[index] ;
	  index <<= 6 ;
	  index += (s[0] << 4) + (s[1] << 2) + s[2] ;
	  ++hex[index] ;
	  ++nHex ;
	}
    }
  fclose (fil) ;

  information (3, codon) ;
  information (6, hex) ;

  if (file2)
    { for (i = 0 ; i < 4096 ; ++i)
	hex2[i] = 1 ;
      nHex2 = 4096 ;

      while (readSequence (fil, dna2indexConv, &seq, &id, 0, &len))
	{ s = seq ;
	  index = (s[0] << 4) + (s[1] << 2) + s[2] ; s += 3 ;
	  for (i = 3 ; i < len-3 ; i += 3, s += 3)
	    { index &= 0x3f ;
	      index <<= 6 ;
	      index += (s[0] << 4) + (s[1] << 2) + s[2] ;
	      ++hex2[index] ;
	      ++nHex2 ;
	    }
	}
      fclose (fil) ;
      hexLikelihoodRatio () ;
    }
  else
    hexTableComposition () ;

  if (ofile)
    saveTable (ofile) ;

  if (sfile)
    { if (!(fil = fopen (sfile, "r")))
	die ("Failed to open fasta file %s", sfile) ;
      nseq = 0 ;
      while (readSequence (fil, dna2indexConv, &seq, &id, 0, &len))
	{ seqs[nseq] = seq ; ids[nseq] = id ; lens[nseq] = len ; 
	  if (++nseq > 1000)
	    die ("More than 1000 sequences - edit and recompile") ;
	}
    }

  scoreSeqs () ;
  return 0 ;
}

/**************** end of file ***************/
