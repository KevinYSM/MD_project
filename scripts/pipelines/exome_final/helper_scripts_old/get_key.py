import sys
import os
import csv

def get_key( var1, var2):
        num_letters=len(var1)
        counter=0
        for i in range(num_letters):
                if (var1[i]!=var2[i]):
                        break
                counter=counter+1
        return var1[:counter]

                

def main():
        key=get_key(sys.argv[1],sys.argv[2])
        print(key)
        return key
main()