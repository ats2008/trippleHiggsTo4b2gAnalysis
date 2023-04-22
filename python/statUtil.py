import numpy as np

def wei_cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - np.average(x, weights=w)) * (y - np.average(y, weights=w))) / np.sum(w)

def wei_corr(x, y, w):
    """Weighted Correlation"""
    return wei_cov(x, y, w) / np.sqrt(wei_cov(x, x, w) * wei_cov(y, y, w))


def twoSample_KSTest(x1,x2,w1=[],w2=[]):
    if w1==[]:
        w1=np.ones(len(x1))
    if w2==[]:
        w2=np.ones(len(x2))

    x1_ = x1 ; w1_ = w1
    x2_ = x2 ; w2_ = w2
    if len(x1) > len(x2):
        x1_=x2 ; w1_=w2
        x2_=x1 ; w2_=w1

    W1=np.sum(w1_)
    W2=np.sum(w2_)
    Dnm=[]
    lb=min(min(x1_),min(x2_))
    ub=max(max(x1_),max(x2_))
    quants=np.linspace( 0.0 , 1.0 , 1000 )
    xq1=np.quantile(x1_,quants)
    xq2=np.quantile(x2_,quants)
    xq=np.concatenate([xq1,xq2])
    for x in xq:
        m1=  x1_  > x 
        m2=  x2_  > x 
        Dnm.append( abs( np.sum(w1_[m1])/W1 - np.sum(w2_[m2])/W2  )  )

    return max(Dnm)

if __name__=='__main__':
    
    x1= np.random.normal(loc=2,scale=3,size=200)
    x2= np.random.normal(loc=3,scale=8,size=200)
    
    rslt=twoSample_KSTest(x1,x2)
    print(rslt)
