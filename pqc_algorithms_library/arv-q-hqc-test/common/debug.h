#ifndef     	__DEBUG_H
#define     	__DEBUG_H
#define         LOG_NONE                0       // no log
#define         LOG_ERROR               1       // printf log level <= LOG_ERROR
#define         LOG_WARNING             2       // printf log level <= LOG_WARNING
#define         LOG_INFO                3       // printf log level <= LOG_INFO
#define         LOG_VERBOSE             4       // printf log level <= LOG_VERBOSE

#define         LOG_LEVEL               LOG_INFO     // set Log level

#define         __BR                      "\r\n"		// line breaker

#if			 	(LOG_LEVEL>=LOG_ERROR)
#define         log_e(format, ...)      printf("E:"format __BR, ##__VA_ARGS__)     
#else
#define         log_e(format, ...)
#endif
#if 			(LOG_LEVEL>=LOG_WARNING)
#define         log_w(format, ...)      printf("W:"format __BR, ##__VA_ARGS__)     
#else
#define         log_w(format, ...)
#endif
#if 			(LOG_LEVEL>=LOG_INFO)
#define         log_i(format, ...)      printf("I:"format __BR, ##__VA_ARGS__)     
#else
#define         log_i(format, ...)
#endif
#if	 			(LOG_LEVEL>=LOG_VERBOSE)
#define         log_v(format, ...)      printf("V:"format __BR, ##__VA_ARGS__)     
#else
#define         log_v(format, ...)
#endif
#endif
