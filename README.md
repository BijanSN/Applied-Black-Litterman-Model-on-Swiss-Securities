# Applied Black-Litterman model on Swiss Equities
A **Quantitative Risk Management** Project

* **GetData.R** : Extract some Swiss bluechips's daily stock prices from yahoo finance and create the datasets.
* **BlackLittermanOnSwissPortfolio.m** : Following predefined specific views on the Swiss market, the associated Optimal Black-Litterman allocation is computed as well as several figures. 

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


# Black-Litterman model

The Black-Litterman method could be interpreted as an extension of the Modern Portfolio
Theory (MPT) since some of the problems that investors have encountered in applying the MPT
(such as high input-sensitivity and mean-variance maximization) are resolved through this
approach. **The idea behind the Black-Litterman model consists of mixing markets assumptions
as well as a manager’s particular view on the market in order to generate a new vector of
expected returns. These views will affect the expected return of the asset concerned, tilting out
our entire portfolio from the efficient MPT frontier.**

The Black litterman model can be divided into 3 differents steps : 

* Get the distribution of the prior
* Get the distribution of the views
* Combining prior and view distributions

## Distribution of the prior

First of all, we combine the regular CAPM's weights of the assets using their respective market capitalisation with the covariance matrix and the risk aversion coefficient in order to compute the implied return vector:

<img src="https://render.githubusercontent.com/render/math?math=\Pi = \lambda  \Sigma W_{Mkt}">

Such as,
* **П**: The implied equilibrium return vector over the risk-free rate
* **λ**: The risk aversion coefficient 
* **Σ**: The Covariance matrix of the assets
* **τ**: The respective weights of views

τ enables more flexibility in our views by introducing uncertainty.
The resulting Prior distribution is therfore normally distributed around П, with variance τΣ. Intuitively, by fixing τ=1, we are fully confident about our views.

## Distribution of the views
The Black-Litterman model uses multiple inputs such as views of a portfolio manager,
expected returns and uncertainty level in order to achieve an optimal asset allocation.

Inputs:

* **P**: The *Pick* matrix that connects the assets related to a specific view.
* **Q**: The vector which represents the expected excess return of each view.
* **Ω**: The diagonal matrix that identifies the uncertainty in the views.

Simply put, the Black-Litterman model is a complex of weighted average of the implied
equilibrium vector (П) and the vector expected excess returns of each view (Q), in which the
relative weights are the scaling factor (τ) and the uncertainty of the forecasts (Ω). Investors
(managers) might have absolute or relative views. In case of absolute view absolute the linked
row in the P-matrix sum up to 1 and in the other case it will sum up to 0. The uncertainty matrix
contains the variance of the views, where ωk is the uncertainty in the k-th view. If ω = 0 ,
the manager if fully confident about the view.

## Combining prior and view distributions

The main formula for implementing this approach is the following :

<img src="https://render.githubusercontent.com/render/math?math=E(R)_{BL}= [ (\tau  \Sigma) ^{-1}  + P^{T}  \Omega^{-1} P]^{-1}  [(\tau   \Sigma )^{-1} \Pi +P^{T} \Omega^{-1} Q]">

Where E[R]BL is the new combined return vector. 

The following illustration is a summary of the mathematical derivation above.
 ![alt text](https://github.com/BijanSN/Applied-Black-Litterman-on-Swiss-portfolios/blob/master/Plots/BL_Summary.PNG)
                                                                                                                 

# Black-litterman model in practice on Swiss equities

WIP

![alt text](https://github.com/BijanSN/Applied-Black-Litterman-on-Swiss-portfolios/blob/master/Plots/PieAllocation.png)

# References

https://ww2.mathworks.cn/help/finance/examples/black-litterman-portfolio-optimization.html
