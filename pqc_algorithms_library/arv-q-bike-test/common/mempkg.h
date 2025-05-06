#ifndef _MEMPKG_H_
#define _MEMPKG_H_
 
#include <stddef.h> /* Using type size_t. */
 
typedef unsigned char UCHART, * PUCHAR;
 
typedef struct st_BLOCK_HEADER
{
    struct st_BLOCK_HEADER * pnext;
    PUCHAR pblock;
    size_t size;
} BLOCK_HEADER, * P_BLOCK_HEADER;
 
#define MEM_SIZ (32768) /* Alter this macro to control the size of a heap. */
 
void   mpkInitMemory(void);
void * Mymemcpy       (void *       dst, void *       src, size_t size);
void * Mymemset       (void *       s,   int          c,   size_t n);
void * Mymalloc       (size_t       size);
void   Myfree         (void *       ptr);
void * Myrealloc      (void *       ptr, size_t       size);
void * Mycalloc       (size_t       n,   size_t       size);
int    Mymemcmp       (const void * s1,  const void * s2,  size_t n);
 
#define memmove(dst, src, size) Mymemcpy(dst, src, size)
 
#endif