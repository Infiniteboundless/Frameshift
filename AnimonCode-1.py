#comparision
def translate(addressIn,randomAmino=[]):
    with open(addressIn,'r')as f:
        seq=f.readline().strip()       
    n=int(len(seq)/3)
    Amino=[x for x in range(0,n)]
    j=0
    for i in range(0,n):
        code=seq[j:3+j]
        j=j+3
        if((code=='UUU')|(code=='UUC')):
            Amino[i]=randomAmino[0]
        elif((code=='UUA')|(code=='UUG')|(code=='CUU')|(code=='CUC')|(code=='CUA')|(code=='CUG')):
            Amino[i]=randomAmino[1]
        elif((code=='AUU')|(code=='AUC')|(code=='AUA')):
            Amino[i]=randomAmino[2]
        elif(code=='AUG'):
            Amino[i]=randomAmino[3]
        elif((code=='GUU')|(code=='GUC')|(code=='GUA')|(code=='GUG')):
            Amino[i]=randomAmino[4]
        elif((code=='UCU')|(code=='UCC')|(code=='UCA')|(code=='UCG')|(code=='AGU')|(code=='AGC')):
            Amino[i]=randomAmino[5]
        elif((code=='CCU')|(code=='CCC')|(code=='CCA')|(code=='CCG')):
            Amino[i]=randomAmino[6]
        elif((code=='ACU')|(code=='ACC')|(code=='ACA')|(code=='ACG')):
            Amino[i]=randomAmino[7]
        elif((code=='GCU')|(code=='GCC')|(code=='GCA')|(code=='GCG')):
            Amino[i]=randomAmino[8]
        elif((code=='UAU')|(code=='UAC')):
            Amino[i]=randomAmino[9]
        elif((code=='UAA')|(code=='UAC')|(code=='UGA')):
            Amino[i]=randomAmino[20]
        elif((code=='CAU')|(code=='CAC')):
            Amino[i]=randomAmino[10]
        elif((code=='CAA')|(code=='CAG')):
            Amino[i]=randomAmino[11]
        elif((code=='AAU')|(code=='AAC')):
            Amino[i]=randomAmino[12]
        elif((code=='AAA')|(code=='AAG')):
            Amino[i]=randomAmino[13]
        elif((code=='GAU')|(code=='GAC')):
            Amino[i]=randomAmino[14]
        elif((code=='GAA')|(code=='GAG')):
            Amino[i]=randomAmino[15]
        elif((code=='UGU')|(code=='UGC')):
            Amino[i]=randomAmino[16]
        elif(code=='UGG'):
            Amino[i]=randomAmino[17]
        elif((code=='CGU')|(code=='CGC')|(code=='CGA')|(code=='CGG')|(code=='AGA')|(code=='AGG')):
            Amino[i]=randomAmino[18]
        elif((code=='GGU')|(code=='GGC')|(code=='GGA')|(code=='GGG')):
            Amino[i]=randomAmino[19]    
    return Amino

#product table
def table(addressIn,addressOut):
    with open(addressOut,'w') as file:
        Aminon=[x for x in'ARDCQEHIGNLKMFPSWTYV*']
        randomAmino=[x for x in range(0,21)]
        flag=0
        for i in range(1,381):
            if(flag>0):
                for j in range(0,20):
                    t=Aminon[0]
                    Aminon[0]=Aminon[j]
                    Aminon[j]=t
                    for n in range(1,20):
                        for n in range(1,20):
                            t=Aminon[n-1]
                            Aminon[n-1]=Aminon[n]
                            Aminon[n]=t
            for j in range(0,21):
                randomAmino[j]=Aminon[j]
            flag=1
            file.write('<')
            for y in translate(addressIn,randomAmino):
                file.write(y)
            file.write('\n')

    
print ('please write down the address of input:')
#addressIn=input()
addressIn="./cytb_cow.fasta" 
print('please write down the address of output:')
#addressOut=input()
addressOut="./sta.txt"
table(addressIn,addressOut)
print ("working...")
