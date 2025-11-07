#include "stdafx.h"

#include <cstring>
#include "formWriter.hpp"

void FormWriter::h1(char *s)
{
	int i,len,L;
	hr();
	hr();
	len = (int)strlen(s)*2;
	L = 40+(frbuf->getIndentLevel()*frbuf->getTabSize()-len)/2;
	for(i=0;i<L;i++) put(' ');
	for(i=0;i<len;i+=2) {put(s[i/2]); put(' ');}
	put('\n');
	hr();
	hr();
}

void FormWriter::h2(char *s)
{
	int i,len,L;
	hr();
	len = (int)strlen(s)*2;
	put('\n');
	L = 40+(frbuf->getIndentLevel()*frbuf->getTabSize()-len)/2;
	for(i=0;i<L;i++) put(' ');
	for(i=0;i<len;i+=2) {put(s[i/2]); put(' ');}
	put('\n');
	put('\n');
	hr();
}

void FormWriter::h3(char *s)
{
	int i,len;
	len = (int)strlen(s);
	for(i=0;i<len;i++) put('-');
	put('\n');
	write(s,len);
	put('\n');
	for(i=0;i<len;i++) put('-');
	put('\n');
}

void FormWriter::center(char *s)
{
	int len,L;
	len = (int)strlen(s);
	L = 40+(frbuf->getIndentLevel()*frbuf->getTabSize()-len)/2;
	write("                                                                               ",L);
	write(s,len);
	put('\n');
}

void FormWriter::hr()
{
	if(!frbuf->preIndent()) {put('\n');}
	write("------------------------------------------------------------------------------",79-frbuf->getIndentLevel()*frbuf->getTabSize());
	put('\n');
	return;
}

int formStreamBuf::sync()
{
   int n = pptr() - pbase();
   return (n && spit(pbase(), n) != n) ? EOF: 0;
}

int formStreamBuf::overflow(int ch)
{
   int n = pptr() - pbase();
   if(n && sync()) return EOF;
   if(ch != EOF) {
      if( spit(ch) != 1) return EOF;
      if( ch == '\n') preindent = 1;
   }
   pbump(-n); // Reset pptr().
   return 0;
}

int formStreamBuf::spit(char *s, int n)
{
   int i,ip=0;
   for(i=0;i<n;i++) {
      if(s[i]=='\n') {
         if(preindent) _indent();
         o.write(&s[ip],i-ip+1);
         ip = i+1;
         preindent = 1;
      }
   }
   if(preindent) _indent();
   o.write(&s[ip],n-ip);
   return n;
}

int formStreamBuf::spit(char c) {
   if(preindent) _indent();
   o.put(c);
   return 1;
}

