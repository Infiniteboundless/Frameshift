#comparision
def translate(addressIn,randomAnimo=[]):
    with open(addressIn,'r')as f:
        seq=f.readline().strip()       
    n=int(len(seq)/3)
    animo=[x for x in range(0,n)]
    j=0
    for i in range(0,n):
        code=seq[j:3+j]
        j=j+3
        if((code=='UUU')|(code=='UUC')):
            animo[i]=randomAnimo[0]
        elif((code=='UUA')|(code=='UUG')|(code=='CUU')|(code=='CUC')|(code=='CUA')|(code=='CUG')):
            animo[i]=randomAnimo[1]
        elif((code=='AUU')|(code=='AUC')|(code=='AUA')):
            animo[i]=randomAnimo[2]
        elif(code=='AUG'):
            animo[i]=randomAnimo[3]
        elif((code=='GUU')|(code=='GUC')|(code=='GUA')|(code=='GUG')):
            animo[i]=randomAnimo[4]
        elif((code=='UCU')|(code=='UCC')|(code=='UCA')|(code=='UCG')|(code=='AGU')|(code=='AGC')):
            animo[i]=randomAnimo[5]
        elif((code=='CCU')|(code=='CCC')|(code=='CCA')|(code=='CCG')):
            animo[i]=randomAnimo[6]
        elif((code=='ACU')|(code=='ACC')|(code=='ACA')|(code=='ACG')):
            animo[i]=randomAnimo[7]
        elif((code=='GCU')|(code=='GCC')|(code=='GCA')|(code=='GCG')):
            animo[i]=randomAnimo[8]
        elif((code=='UAU')|(code=='UAC')):
            animo[i]=randomAnimo[9]
        elif((code=='UAA')|(code=='UAC')|(code=='UGA')):
            animo[i]=randomAnimo[20]
        elif((code=='CAU')|(code=='CAC')):
            animo[i]=randomAnimo[10]
        elif((code=='CAA')|(code=='CAG')):
            animo[i]=randomAnimo[11]
        elif((code=='AAU')|(code=='AAC')):
            animo[i]=randomAnimo[12]
        elif((code=='AAA')|(code=='AAG')):
            animo[i]=randomAnimo[13]
        elif((code=='GAU')|(code=='GAC')):
            animo[i]=randomAnimo[14]
        elif((code=='GAA')|(code=='GAG')):
            animo[i]=randomAnimo[15]
        elif((code=='UGU')|(code=='UGC')):
            animo[i]=randomAnimo[16]
        elif(code=='UGG'):
            animo[i]=randomAnimo[17]
        elif((code=='CGU')|(code=='CGC')|(code=='CGA')|(code=='CGG')|(code=='AGA')|(code=='AGG')):
            animo[i]=randomAnimo[18]
        elif((code=='GGU')|(code=='GGC')|(code=='GGA')|(code=='GGG')):
            animo[i]=randomAnimo[19]    
    return animo

#product table
def table(addressIn,addressOut):
    with open(addressOut,'w') as file:
        animon=[x for x in'ARDCQEHIGNLKMFPSWTYV*']
        randomAnimo=[x for x in range(0,21)]
        flag=0
        for i in range(1,381):
            if(flag>0):
                for j in range(0,20):
                    t=animon[0]
                    animon[0]=animon[j]
                    animon[j]=t
                    for n in range(1,20):
                        for n in range(1,20):
                            t=animon[n-1]
                            animon[n-1]=animon[n]
                            animon[n]=t
            for j in range(0,21):
                randomAnimo[j]=animon[j]
            flag=1
            file.write('<')
            x=translate(addressIn,randomAnimo)
            for y in x:
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
