#!/usr/bin/env python
#
#  word_wrap.py
#  
#
#  Created by Danny Jacobs on 6/10/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#
import numpy as n

def word_wrap(string, width=80,ind1=0,ind2=0,prefix=''):
    """ word wrapping function.
        string: the string to wrap
        width: the column number to wrap at
        prefix: prefix each line with this string (goes before any indentation)
        ind1: number of characters to indent the first line
        ind2: number of characters to indent the rest of the lines
    """
    awidth = min(width-2-len(prefix+ind1*' '),width-2-len(prefix+ind2*' '))
    outstring = ''
    words = string.split(' ')
    okwords = []
    chunk = lambda v, l: [v[i*l:(i+1)*l] for i in range(int(n.ceil(len(v)/float(l))))]
    for word in words:
        for okword in chunk(word,awidth):
            okwords.append(okword)
    lines = []
    l = prefix+ind1*' '
    for i,w in enumerate(okwords):
        print w,len(l+' '+w),width
        if len(l+' ' + w)<width:
            l += ' '+w
        else:
            lines.append(l)
            l = prefix + ind2*' '+w
    lines.append(l)
            
    return '\n'.join(lines)+'\n'
if __name__=='__main__':
    teststr1 = 'stupidlongfilename'*5 + ' other stuff '
    print teststr1
    print '-'*10
    print word_wrap(teststr1,70,5,20,'#')     
            
    