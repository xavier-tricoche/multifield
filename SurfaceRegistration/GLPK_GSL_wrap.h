#include <gsl/gsl_vector.h>	/*GNU Scientific Library*/
#include <gsl/gsl_matrix.h>
#include <glpk.h>		/*GNU Linear Programming Kit*/

LPX *setUpLOneNormMinimisation(gsl_matrix *);
void LOneNormMinimisationVector(LPX *lp, gsl_vector *, gsl_vector *);
void LOneNormMinimisationMatrix(gsl_matrix *, gsl_matrix *, gsl_matrix *);
