import numpy as np 
import scipy.stats as stat
import matplotlib.pyplot as plt

data=np.load("PTZ-WILDTYPE-02_2photon_sess-01-6dpf_BLN_run-01_0.590bin0.10nnbav.npy")
sizes=data[0,:]
M=len(sizes)
a=min(sizes)
b=max(sizes)

def powerlaw(n,lam):
    zeta=np.sum(1.0/np.arange(a,b+1)**lam)
    return(n**(-lam)/zeta)

def lognormal(n,mu,sig):
    return(1.0/n/np.sqrt(2*np.pi*sig**2)*np.exp(-(np.log(n)-mu)**2/(2*sig**2)))

def LogLikelihood(lam):
    zetamat=np.power.outer(1.0/np.arange(a,b+1),lam)
    zeta=np.sum(zetamat,0)
    nprod=-lam*np.sum(np.log(sizes))
    norm=-M*np.log(zeta)
    loglik=nprod+norm
    return(loglik) 

def LogLikelihood_LN(mu,sig):
    T1 = -np.sum(np.log(sizes))
    T2_mat = np.subtract.outer(np.log(sizes),mu)**2
    T2 = -np.sum(T2_mat,0)/(2*sig**2)
    T0 = -M*np.log(np.sqrt(2*np.pi) * sig )
    loglik=T0+T1+T2
    return(loglik) 

def IS(npart):
    lambda_sample=np.random.uniform(0.1,5,npart)
    weights=LogLikelihood(lambda_sample)+stat.norm.logpdf(lambda_sample,1,3)-stat.uniform.logpdf(lambda_sample,0.1,5)
    maxw=np.max(weights)

    w2 = np.exp(weights-maxw)
    w2_sum = np.sum(w2)

    ESS=1.0/(np.sum((w2/w2_sum)**2))

    mean_lambda = np.dot(lambda_sample,w2)/w2_sum
    return([mean_lambda,maxw + np.log(np.sum(np.exp(weights-maxw)))-np.log(npart),ESS])

def IS_LN(npart):
    mu_sample = np.random.uniform(-2.0,2.0,npart)
    sig_sample = np.random.uniform(0.1,5.0,npart)
    weights=LogLikelihood_LN(mu_sample,sig_sample)
    maxw=np.max(weights)

    w2 = np.exp(weights-maxw)
    w2_sum = np.sum(w2)

    ESS=1.0/(np.sum((w2/w2_sum)**2))

    wmax_ID=np.argmax(w2)
    mean_mu = mu_sample[wmax_ID]
    mean_sig = sig_sample[wmax_ID]

    return([mean_mu,mean_sig, maxw + np.log(np.sum(np.exp(weights-maxw))),ESS])


def plot_samples(npart):
    lambda_sample=np.random.uniform(0.1,5,npart)
    weights=LogLikelihood(lambda_sample)
    maxw=np.max(weights)
    w2 = np.exp(weights-maxw)
    
    plt.hist(lambda_sample,weights=w2,bins=np.linspace(2.5,2.8))
    plt.show()

def plotcomp(lam,mu,sig):
    x = np.linspace(a,b,40) 
    plt.hist(sizes,40,log=True,density=True)
    plt.plot(x,powerlaw(x,lam))
    plt.plot(x,lognormal(x,mu,sig))
    plt.show()

ln=IS_LN(2000)
po=IS(2000)
plotcomp(po[0],ln[0],ln[1])
