#ifndef STDLIBS
#define STDLIBS                                                                                                                                                                              
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <stdarg.h>
#include <sys/time.h>
#endif

#define NTAXA         400      /* max # of species */
#define NGENE         20000      /* max # of loci */
#define MAXROUND	10000000		/* MAX # OF ROUNDS*/
#define NUM_NOCHANGE	20000		/* # OF ROUNDS THAT NO BIGGER LIKELIHOOD VALUES ARE FOUND*/
#define LSPNAME       60       /* # characters in sequence names */
#define MAXBRLENS 7.0 /* the maximum branch length of the species tree */
#define COLLAPSEBRLENS 1e-6 /* the branch length for collapsing branches */  
#define ERROR 1
#define NO_ERROR 0
#define YES 1
#define NO 0
#define NA -1
#define DEBUG 0
#define FPN(file) fputc('\n', file)
#define FOR(i,n) for(i=0; i<n; i++)
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))

typedef struct node 
	{
	int father, nson, sons[2], namenumber;
	char taxaname[LSPNAME];
	double brlens,theta;
   	}
	Treenode;

typedef struct Tree
	{
   	int root;
	int ntaxa; 
   	Treenode nodes[2*NTAXA];
	}  
	Tree;

/* tool functions*/
FILE *gfopen(char *filename, char *mode);
void SetSeed (unsigned int seed);
double LnGamma (double x);
double rndu (void);

int 		addNode (Tree *tree, int fromnode, int fromfather, int tonode);
int 		AddToPrintString (char *tempStr);
int 		Algorithm (int **triple, FILE *outfile, FILE *outputfile);
void 		copyTree(Tree *from, Tree *to);
int 		deleteNode (Tree *tree, int inode);
int 		genetreeTriples (Tree *tree, FILE *outfile);
int 		findNgenesandNtaxa (FILE *fTree);
int 		findNameNumber (Tree *tree);
int 		findTriple (int node1, int node2, int node3, long int *location);
void 		findOffsprings (int *offsprings, Tree *tree, int inode);
int 		findOutgroup (Tree *tree);
int 		logbinomialP (int n, int x, double p, double *logp);
int 		logLikelihood (int **triple, Tree *tree, double *loglike);
int 		maximizeaBrlen (int **triple, Tree *tree);
void 		MoveBrlens (Tree *tree);
int 		MoveNode (Tree *tree);
void		MrBayesPrint (char *format, ...);
void 		PrintHeader (void);
int 		printSptree (void);
int 		PrintState (int round, FILE *fout, int addend);
int 		PrintTree (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void 		printTriples (FILE *outfile);
void 		randomTree(Tree *tree);
void 		randomVector (int *array, int number);
int 		ReadaTree (FILE *fTree,Tree *tree);
void 		*SafeMalloc(size_t s);
int 		SaveSprintf(char **target, int *targetLen, char *fmt, ...);
int 		swapNodes (Tree *tree, int inode, int jnode);
int 		tripleDist (Tree *tree1, Tree *tree2);
int 		triples (Tree *tree, int ntrees);
int 		triplesInatree (int **triple, Tree *tree);
int			triplesInatree1	(int **triple, Tree *tree);
void 		WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void		WritePhylipTreeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);
int			PrintPhylipTree(Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);



