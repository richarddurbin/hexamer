/*  File: readseq.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, Wellcome Sanger Institute, 1993-1998
 *  License: available under the MIT license as in the accompanying LICENSE file 
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  2 23:50 2021 (rd109)
 * Created: Tue Jan 19 21:14:35 1993 (rd)
 *-------------------------------------------------------------------
 */

extern int readSequence (FILE *fil, int *conv,
			 char **seq, char **id, char **desc, int *length) ;
				/* read next sequence from file */
extern int writeSequence (FILE *fil, int *conv, 
			  char *seq, char *id, char *desc, int len) ;
				/* write sequence to file, using convert */
extern int seqConvert (char *seq, int *length, int *conv) ;
				/* convert in place - can shorten */
extern int readMatrix (char *name, int *conv, int** *mat) ;

extern int dna2textConv[] ;
extern int dna2indexConv[] ;
extern int dna2binaryConv[] ;
static const char index2char[] = "acgtn" ;
extern int aa2textConv[] ;
extern int aa2indexConv[] ;
static const char index2aa[] = "ACDEFGHIKLMNPQRSTVWYX" ;

/***** end of file *****/
