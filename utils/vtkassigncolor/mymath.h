#if defined(WIN32) || defined(__WIN32) || defined(__WIN32__) || defined(WIN64) || defined(__WIN64) || defined(__WIN64__) 
#define ISNAN _isnan
#else
#define ISNAN isnan
#endif

#define MFL { printf("%d in %s\n", __LINE__, __FILE__); }