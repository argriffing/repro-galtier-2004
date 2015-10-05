from __future__ import print_function, division

import argparse


templated_options = """\
/* all programs */

  /*
   * 1 is optimized, 0 is not optimized
   */
  OPTIMIZE_LENGTH=1		/* branch lengths */
  OPTIMIZE_GC=0			/* G+C contents */
  OPTIMIZE_TITV=1		/* transition/transversion ratio */
  OPTIMIZE_ROOT=1		/* root location */
  OPTIMIZE_ANC=0		/* ancestral G+C content */
  OPTIMIZE_GAMMA=1		/* shape of the gamma distribution */
  OPTIMIZE_COV=1		/* covarion rate */
  OPTIMIZE_PI=0			/* proportion of covarion sites */

  INIT_LENGTH=REDO 		/* KEEP or REDO */
  INIT_GC=CONST			/* CONST or VAR or BALANCED */
  INIT_TITV=4.			/* -1.=automatic */
  INIT_ROOT=0.5 		/* -1.=keep */
  INIT_ANC=-1.                  /* -1.=automatic */
  INIT_GAMMA=0.3		/* -1.=infinity (constant rate) */
  GAMMA_NBCLASS={GAMMA_NBCLASS}		/* >= 2 */	
  INIT_COV=0.1			/* -1.=infinity (no covarion) */
  INIT_PI=0.			
		
  PRECISION=4
  SIMPLEX=0

  PRINT1=1
  PRINT2=1



/* eval */

  ALLCOUPLES=1


/* shake*/

  SH_G=-1				/* -1=max */
  SH_MAXLCROSSED=-1.			/* -1.=> every branch crossed */
  SH_MAXBOOTCROSSED=90			/* -1=> every branch crossed */
"""


def nbclass_type(s):
    g = int(s)
    if g < 2:
        raise ValueError('each number of gamma classes must be at least 2')
    return g


def main(args):
    for g in args.nbclass:
        filename = '{}{}{}'.format(args.prefix, g, args.suffix)
        with open(filename, 'w') as fout:
            contents = templated_options.format(GAMMA_NBCLASS=g)
            print(contents, file=fout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', default='nhml-options-')
    parser.add_argument('--suffix', default='.opt')
    parser.add_argument('nbclass', type=nbclass_type, nargs='+',
            help='GAMMA_NBCLASS to be passed to nhml')
    args = parser.parse_args()
    main(args)

