# Applied Black-Litterman model and Copula Opinion Pooling on Swiss Equities
A **Quantitative Risk Management** Project

# Introduction

Departing from the CAPM’s strong assumption of the homogeneity of expectations of risk
and returns, particular views of a manager can be implemented into the asset allocation, using
the Black-Litterman model (BL) in order to potentially generate alpha, the sought outperformance over
the benchmark. After a brief review of this model and its assumptions, a basic implementation
on Swiss equities will be given.
Then, the COP (the copula-opinion pooling) will be presented and a succinct
example on the yield curve will be provided. **The underlying question of the project will be to
know if our particular views indeed help to achieve a better allocation under this particular model**. In
conclusion, we will review the results, discuss new trends in finance and emphasize why using
the above-mentioned approach is relevant nowadays.


# Modern Portfolio Theory deficiencies
Through a formalisation of the benefits of diversification and a constrained optimisation using two-dimensions (risk and returns), the Modern Portfolio Theory (MPT) - developed by Harry Markowitz- demonstrated that among the wide spectrum of securities, an optimal combinaison of securities can be found for each level of risk or returns.
As such, the expected return can be maximized for a given level of risk and vice versa. Those combinaisons, located on the efficient frontier, are called efficient portfolios.
Among the infinite combinaisons of portfolios on the efficient frontier, some notorious ones are the **Minimum Variance (MVP)** and the **Maximum Sharpe Ratio (MSR)** portfolios. The MVP is simply the portfolio with the lowest risk, the first on the efficient frontier. The MSR portfolio is the portfolio which it's Sharpe ratio (expected excess returns over it's variance) is the highest among the efficient frontier.

However, this framework receives many criticisms, such as simplifying assumptions (the Gaussian distribution of the returns for example), the lack of consideration for risk asymetry during downturns and the high input sensitivity, leading very different results given a small change in inputs. The latter problem will be dealt with the COP approach onward. But most importantly,**the MPT considers that each investors share similar views on expected returns on securities**. One could argue that each manager have valuable knowledge on their specific area of expertise, leading , hopefully, in the creation of alpha. As an introduction, the Black litterman model presented below is mearly a weighted average of these valuables views and the markets inputs.

To overcome the MPT deficiencies and provide better flexibility for the manager, Fisher Black and Robert Litterman developed in 1990 the Black Litterman (BL) asset allocation model. The BL approach allows to combine the MPT with some subjectivity, intuitions, particular views of the investor to the equation. It is important to stress that this newly-created portfolio is, following the traditional MPT, below the efficient frontier and therefore considered irrational because we could gain more returns for less risks by choosing the efficient one given the same level of risk.


## 2.1 Black-Litterman with a normally distributed prior

The Black-Litterman method could be interpreted as an extension of the Modern Portfolio
Theory (MPT) since some of the problems that investors have encountered in applying the MPT
(such as high input-sensitivity and mean-variance maximization) are resolved through this
approach. **The idea behind the Black-Litterman model consists of mixing markets assumptions
as well as a manager’s particular view on the market in order to generate a new vector of
expected returns. These views will affect the expected return of the asset concern, tilting out
our entire portfolio from the efficient MPT frontier.**

The Black litterman method can be divided into 3 differents steps : 

* Get the distribution of the prior
* Get the distribution of the views
* Combining prior and view distributions


# Distribution of the prior
First of all, we combine the regular CAPM's weights of the assets using their respective market capitalisation with the covariance matrix and the risk aversion coefficient in order to compute the implied return vector:

<img src="https://render.githubusercontent.com/render/math?math=\Pi = \lambda  \Sigma W_{Mkt}">

Such as,
* П The implied equilibrium return vector over the risk-free rate
* λ The risk aversion coefficient 
* Σ The Covariance matrix of the assets
* τ The respective weights of views

τ enables more flexibility in our views by introducing uncertainty.
The resulting Prior distribution is therfore normally distributed around П, with variance τΣ. Intuitively, by fixing τ=1, we are fully confident about our views.


### 2. Distribution of the views
The Black-Litterman model uses multiple inputs such as views of a portfolio manager,
expected returns and uncertainty level in order to achieve an optimal asset allocation.

Inputs:

*P The “Pick” matrix that connects the assets related to a specific view.
*Q The vector which represents the expected excess return of each view.
*Ω The diagonal matrix that identifies the uncertainty in the views.


Simply put, the Black-Litterman model is a complex of weighted average of the implied
equilibrium vector (П) and the vector expected excess returns of each view (Q), in which the
relative weights are the scaling factor (τ) and the uncertainty of the forecasts (Ω). Investors
(managers) might have absolute or relative views. In case of absolute view absolute the linked
row in the P-matrix sum up to 1 and in the other case it will sum up to 0. The uncertainty matrix
contains the variance of the views, where ωk is the uncertainty in the k-th view. If ω = 0 ,
the manager if fully confident about the view.


### 3. Combining prior and view distributions

The main formula for implementing this approach is the following :
<img src="https://render.githubusercontent.com/render/math?math=E(R)_{BL}= [  (\tau  \Sigma) ^{-1}  + P^{T}  \Omega^{-1} P]^{-1}  [(\tau   \Sigma )^{-1} \Pi +P^{T} \Omega^{-1} Q]">


<img src="https://render.githubusercontent.com/render/math?math=E(R)= [  (\tau  \Sigma) ^{-1}  + [P^{T}  \Omega P]^{-1}  [(\tau   \Sigma )^{-1} \Pi +P^{T} \Omega Q]">

Where,
E[R]BL is the new combined return vector 


The following illustration is a summary of the mathematical derivation above.


# 1. Black-Litterman using copula approach (COP)

One of the main drawbacks of the previous BL implementation is embedded into the former
prior described. Indeed, as stated earlier, **the prior follows a normal distribution around the
mean of the excessive market returns. In reality, this is not the case at all, as our empirical data
suggests. To bypass this problem, the copula-opinion pooling (COP) - an extension in a classical
Black-Litterman method- implements a prior which correctly reflects the dependency structure
of the market using multivariate copula**. The intuition is the same as before; the Black-Litterman
method requires 2 inputs: Market views and well as the manager’s. In order to infer the nonparametric
market prior this time, the probability density function of the market should be
estimated. The intuition being that past distribution of the market reflects our “prior”
assumption on how the market will behave in the future, which seems as a more reasonable
assumption than a normal distribution around the excess market returns.
The prior market distribution is computed through Monte Carlo simulations, drawn from
the previous prior distribution computed before. Those simulations are not subject to the curse
of dimensionality, and represent an efficient way to construct market views. We can now depart
from the standard efficient allocation by integrating some of our subjective views on the market.
The second input of the Black-Litterman method is the manager’s views, which is left
unchanged in this case. As before, the posterior is defined as a weighted average between the
market and the views, with a confidence coefficient associated to both.

## 3.1 Example on US treasury yield
In this section, we will highlight the drawbacks of using a normally distributed prior using
empirical data. Using data on the US treasury yield, the COP approach has been implemented,
using “Beyond Black-Litterman in Practice:a Five-Step Recipe to Input Views on non-Normal
Markets” from Attilio Meucci, as a benchmark for our analysis. Part of the matlab code of the
project are borrowed from him in this section. The yield rate are plotted below from different
maturities, from 3 months to 30 years:


The plot is coherent with the economic theory; ceteris paribus, for a higher maturity
corresponds to a higher risk premium and therefore a higher yield. Inversions of the yield curve
may suggest worsening expectation concerning the future, hence higher maturities carrying
lesser yield. Furthermore, the yield curve can be derived as well, adding a third dimension:
maturities.

We can extract the linear dependency among the curves by the simple pearson correlation,
as shown with this heatmap:

As our intuition suggests, the closer on yield is from another one, the higher is the
correlation between them, since they are affected by the same risk factors. The correlation still
remains above 0.6, which is consequent. Now, we want to capture other form of dependences,
through copula. As described before, the first step consists of extracting the market views.
The empirical (unsmoothed) CDF is shown below. Indeed, the plot looks as a “stair”
function, which didn’t allow the derivation into the PDF. A kernel smoothing is applied on the
CDF in order to correctly represent the density of the change in yield.

Compute the market distribution now yield the following results:

As we can see, the distribution of the yield changes are leptokurtic and mostly centered
around zero. Applying univariate Jarque-Bera tests and Mardia’s multivariate test on the change
of yield shows a clear rejection from the normality distribution assumption mostly due to the
presence of these fat negative tails. Graphically, looking specifically at the 3 months yield
change and a gaussian distribution shows a clear difference on the distribution of the yield than
the traditional Black-Litterman model assumption, hence the importance of choosing a
different market prior.

We can now draw simulations of those distribution in order to simulate our market prior.
For example, the histogram of the 3 month yield versus the 2 year yield can be shown below.
It is not really informative by itself, but by extracting the margins from theses distribution, the
notion of copula emerges.

As we can see, there is a clustering of data on both side of the tails, which corresponds to
a t-copula. We will therefore use this model to represent the codependency of the weekly
changes. Finally, the parameter of the t-copula are estimated through a maximum-likelihood.
Following the same example as Attilio Meucci (2006), the views are defined such that the
investor expects a steepening of the intermediate yield curve. More specifically, the manager
expected an increase of 10 basis point of the spread 2-20 years as well as a 5 .bp bullish view
on the 2-5-10 butterfly (short 2th&10th year and long the 5th year yield)
20
The first line of the pick matrix is therefore [0 -1 0 0 0 1] with an associated coefficient of
0.001. The second line is [0 -0.5 1 -0.5 0 0] with a coefficient of 0.0005. The vector of expected
excess returns of each view ‘Q’ is therefore [0.001, 0.0005].
Intuitively, higher yield leads to a decrease of the bonds price, the manager should have a
short position on bonds for the 5 year yield as well the 20 year yield, according to his view. For
the second view, the proceed of the sell of the 5 year maturity is used on both the 2 year and 10
year maturity to bet on a steepening of the yield curve around that maturity.
From this plot, we can see that the prior market assumption is blended with our particular
view, in order to create a new, posterior view on the market. This posterior view can be
implemented in order to create a new allocation, depending on the framework we are using.
A summary table of the difference between the prior and the posterior yield change are
given below: Each column represent the different maturities, respectively 2,5,10 and 20 years.

As we can expect,the view of the manager twisted the expected mean returns of the
posterior in line with his expectations. Indeed, the expected mean of the change in yield for the
2nd maturity of the posterior is slightly lower than the prior’s, consequence of the both views
of the manager, the fact that the slope of the yield curve is expected to steepen. Similarly,the
10th yield increases as we expect the spear 2-20 to widen. As for the 5th yield, the bullish view
is represented by a higher expected yield change. For the standard deviations , the skewness,
and the kurtosis, no significant changes arises between the prior and posterior expected yield
change. However, we can easily see that for both the prior and posterior, the overall skewness
of the change in yield is slightly negative and the overall kurtosis is above 3, indicating a
leptokurtic distributions of the yield change, indicating again a departure from the gaussian
approximation.


##2. Conclusion
We are at the edge of a financial revolution. With the surge of data, machine learning and
artificial intelligence techniques applied to finance are becoming more efficient at solving
financial-related problems than us, humans. A legitimate controversial question can be brought
here; Are we still going to be useful? Are we going to be useful in the future? What is the
competitive advantage a human has over an algorithm? We believe it’s our intuition. The
Black-Litterman model combines these two paradigms; By combining the expertise we
gathered the last couples of years with our intuition, we can achieve better results, as seen
through the swiss equity example. The assumption of normality was also rejected several times
through our project, leading the development of new tools in order to correctly model financial
information. This was highlighted through the COP approach from Meucci, but the research
must go on. To conclude with, we presented a relatively new approach to asset management :
combining a specialist view with the market’s can potentially yields higher results, depending
on certains criteria, such as the confidence of the manager on his views. An extension would
be to analyse if an optimal confidence level (yielding a better allocation) could be inferred

#Reference 
[1](http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=5472145&fileOId=5472158)





