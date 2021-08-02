/*  File: hexamer.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) 1993-2021 Richard Durbin, Wellcome Sanger Institute
 *  License: available under the MIT license as in the accompanying LICENSE file 
 *-------------------------------------------------------------------
 * Description: finds maximal scoring coding potential segments
                based on log likelihood ratio scores made by hextable
		of coding hexamers to all the hexamers with the same base composition
 * Exported functions: main()
 * HISTORY:
 * Last edited: Aug  2 23:49 2021 (rd109)
 * * Aug  2 22:59 2021 (rd109): removed all acedb code in this standalone version
 * Created: Sun Aug 27 16:08:28 1995 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>		/* defines bool, true, false */
#include <string.h>

/*-----------------------------------------------------------*/

static bool readTable (char *name, float **tab)
{
  int j, n = 0 ;
  FILE *fil ;
  int level ;

  if (!(fil = fopen (name, "r"))) return false ;

  if (!*tab) *tab = (float*) malloc (4096*sizeof(float)) ;

  for (j = 0 ; j < 4096 && !feof(fil) ; ++j)
    if (fscanf (fil, "%f", &(*tab)[n]))
      ++n ;
    else
      { fprintf (stderr, "can't find entry %d in table file\n", j) ;
        return false ;
      }

  if (n == 4096)
    return true ;
  else
    { fprintf (stderr, "problem reading hexamer file %s.hex", name) ;
      return false ;
    }
}

/********** form partial sum array ***********/

static float makePartial (char *s, int len, float *tab, int skip, float *partial)
{
  int i, j, index ;
  float score = 0 ;

  if (len < 6) return 0. ;

  index = (s[0] << 10) + (s[1] << 8) + (s[2] << 6) + (s[3] << 4) + (s[4] << 2) + s[5] ;
  s += 6 ;
  for (i = 3 ; i <= len-3 ; i += skip)
    { score += tab[index] ;
      partial[i] = score ;
      for (j = skip ; j-- ;) index = (index << 2) + *s++ ;
      index &= 0xfff ;
    }

  return score ;
}

/***** find and print maximal segments *****/

static void printSeg (int x1, int x2, float score) ;

static void processPartial (int step, float thresh, bool isRC,
			    int offset, float *partial, int len)
{ 
  register int i, k ;
  int loclen ;
  static int *maxes = 0, *mins = 0, slen = 0 ;

  loclen = len - offset ;
  while (loclen % step) --loclen ;
  partial += offset ;

  if (slen < len)
    { free (maxes) ; maxes = (int*) malloc (len * sizeof(int)) ;
      free (mins) ; mins = (int*) malloc (len * sizeof(int)) ;
      slen = len ;
    }

  k = 3 ;					/* make mins */
  for (i = 3 ; i <= loclen-3 ; i += step)
    { if (partial[i] < partial[k]) k = i ;
      mins[i] = k ;
    }
  k = loclen - 3 ;				/* make maxes */
  for (i = loclen - 3 ; i >= 3 ; i -= step)
    { if (partial[i] > partial[k]) k = i ;
      maxes[i] = k ;
    }

  for (i = 3 ; i <= loclen - 3 ; i += step)
    if (mins[maxes[i]] == i && 
	partial[maxes[i]] - partial[i] > thresh)
      { if (isRC)
	  printSeg (len-1 - maxes[i] - offset, len-1 - i - offset,
		    partial[maxes[i]] - partial[i]) ;
	else
	  printSeg (i+offset, maxes[i]+offset, partial[maxes[i]] - partial[i]) ;
      }
}

/****************************************************************/

#include "readseq.h"

/****************************************************************
** stand alone version
** writes a gff file to stdout
*****************************************************************/

static void usage (void)
{
  fprintf (stdout, "Usage: hexamer [opts] <tableFile> <seqFile>\n") ;
  fprintf (stdout, "options: -T <threshold>          0\n") ;
  fprintf (stdout, "         -F <feature name>       tableFile name\n") ;
  fprintf (stdout, "         -n	   flag for noncoding (no triplet frame)\n") ;
  exit (-1) ;
}

static char *seqName = "" ;
static char *tableName ;
static char *featName ;
static char strand ;
static char frame = '0' ;

static void printSeg (int x1, int x2, float score)
{
  printf ("%s\t%s\t%s\t%d\t%d\t%.4f\t%c\t%c\n", 
	  seqName, "hexamer", featName, x1+1, x2+1, 
	  score, strand, frame) ;
}

int main (int argc, char **argv)
{
  int i ;
  float thresh = 0.0 ;
  float *tab = 0, *partial ;
  int step = 3 ;
  FILE *seqFile ;
  char *seq ;
  char c ;
  int len ;
  bool isRC ;

  --argc ; ++argv ;		/* remove program name */

  while (argc > 2)
    if (!strcmp (*argv, "-T"))
      { thresh = atof (argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-F"))
      { featName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-n"))
      { step = 1 ; frame = '.' ;
	argc -= 1 ; argv += 1 ;
      }
    else if (**argv == '-')
      { fprintf (stderr, "Unrecognised option %s\n", *argv) ;
	usage() ;
      }
    else
      usage() ;

  if (argc != 2)
    usage() ;

  tableName = *argv ; --argc ; ++argv ;
  if (!featName) featName = tableName ;
  if (!readTable (tableName, &tab))
    { fprintf (stderr, "Failed to open table file %s\n", tableName) ;
      usage () ;
    }

  if (!strcmp (*argv, "-"))
    seqFile = stdin ;
  else if (!(seqFile = fopen (*argv, "r")))
    { fprintf (stderr, "Failed to open sequence file %s\n", *argv) ;
      usage() ;
    }

  if (!readSequence (seqFile, dna2indexConv, 
		     &seq, &seqName, 0, &len))
    { fprintf (stderr, "Errors reading sequence %s\n", *argv) ;
      usage() ;
    }
  for (i = 0 ; i < len ; ++i)	/* remove N's - map to C for this */
    if (seq[i] > 3) seq[i] = 1 ;

  partial = (float*) malloc (sizeof(float)*len) ;

				/* first do forward direction */
  isRC = false ; strand = '+' ;
  for (i = 0 ; i < step ; ++i)
    makePartial (seq+i, len-i, tab, step, partial+i) ;

  for (i = 0 ; i < step ; ++i)
    processPartial (step, thresh, isRC, i, partial, len) ;

				/* then reverse complement */
  isRC = true ; strand = '-' ;
  for (i = 0 ; i < len-1-i ; ++i)
    { c = 3 - seq[i] ;	/* NB "3 -" does complement */
      seq[i] = 3 - seq[len-1-i] ; 
      seq[len-1-i] = c ; 
    }

  for (i = 0 ; i < step ; ++i)
    makePartial (seq+i, len-i, tab, step, partial+i) ;
  for (i = 0 ; i < step ; ++i)
    processPartial (step, thresh, isRC, i, partial, len) ;

  free (seq) ;
  free (partial) ;
}

/**************** end of file ****************/
