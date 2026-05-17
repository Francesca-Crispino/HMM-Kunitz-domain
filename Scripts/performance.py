#!/urs/bin/LAB1 python
import sys
import numpy as np

def get_preds(fname): # "get predictions"
    preds = []
    fh = open(fname)

    for line in fh:
        v = line.rstrip().split()  # extracting the name of the protein, name of the best match and the e-value
        preds.append([v[0], float(v[1]), int(v[2])]) 
    return preds

def get_cm(preds, th = 0.001): # 'get confusion matrix'
    cm = np.zeros((2,2))
    n = len(preds)
    for k in range(n):
        j = 0
        i = preds[k][2]
        
        if len(preds[k]) < 3: continue
        if float(preds[k][1]) <= th: 
            #Compare the sequence E-value against the selected threshold. If the E-value is smaller than or equal to the threshold,
            #the sequence is classified as positive (j = 1), otherwise it is classified as negative (j = 0).
            j=1 
        else: j=0
        cm[i,j]=cm[i,j]+1
    return cm


def get_acc(cm): # 'get accuracy'
    total=np.sum(cm)
    return (cm[0,0] + cm[1,1])/total 


def get_mcc(cm):   #Matthew Correlation Coefficient 
    tp = cm [1,1]  #True Positives
    tn = cm[0,0]   #True Negatives
    fn = cm[1,0]   #False Negatives 
    fp = cm[0,1]   #False Positives 
    numerator = (tp *tn)-(fp*fn) #Calculate the denominator 
    d=(tp + fp) * (tp + fn) * (tn +fp) * (tn +fn)
    mcc=numerator/np.sqrt(d)
    return mcc

if __name__ == '__main__':
    fname = sys.argv[1]   
    th = float(sys.argv[2])      # threshold
    preds = get_preds(fname)
    cm = get_cm(preds,th)
    q2 = get_acc(cm)
    mcc = get_mcc(cm)
    print("TH:",th , "Q2:",q2,"MCC",mcc)
    

    



