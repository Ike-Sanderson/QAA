#!/usr/bin/env python

#version 1.1

maps = 0
unmap = 0

file = "/projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_24_4A_Aligned.out.sam"
with open (file, 'r') as fh:
    for line in fh:
        #print(line)
        if line [0] != '@':
            #print(line)
            flag = int(line.split()[1])
            if((flag & 256) != 256): #check for secondary alignments
                if((flag & 4) != 4): #mapped or unmapped
                    #mapped = True
                    maps += 1
                else:
                    #mapped = False
                    unmap += 1
            #else:                    #secondary reads count as mapped, NOT unmapped
                #mapped = False
                #unmap += 1
#print(maps, unmap)
print(f'{maps} reads were mapped')
print(f'{unmap} reads were unmapped')