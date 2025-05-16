# 🎲 CMPE320 Project 3 – The Central Limit Theorem: The Magic of Large Numbers

## 📌 About

In this project, we explore one of the most important ideas in probability: the Central Limit Theorem (CLT). We’ve heard that when you add up a bunch of random variables—even if they come from different distributions—the result starts to look like a normal (Gaussian) distribution. Here, we actually test that out ourselves.

We simulate thousands of trials where we take random variables, add them together, and plot histograms of the results. We repeat this for different kinds of distributions: uniform, discrete (like dice rolls), exponential-like, and even binary ones. We see how the shape of the resulting histogram changes depending on how many variables we add and what kind of distribution we started with. We also calculate the sample mean and variance to compare with the expected values, and we plot theoretical Gaussian curves to visually compare.

Everything is done in MATLAB, and we use what we learned from Projects 1 and 2 to create histograms, overlay PDFs, and interpret the results.

---

## 🎯 Goal

The goal is to:
- Understand and demonstrate how the Central Limit Theorem works.
- Simulate the sum of independent, identically distributed (i.i.d.) random variables from different distributions.
- Compare the sample results with theoretical expectations for mean and variance.
- Overlay Gaussian curves to show how the distributions approach normality as the number of terms increases.
- Analyze both continuous and discrete random variables to show how CLT applies in each case.

---

## ⚙️ How It Works

Here’s what we do:
1. **We simulate sums of uniform random variables** (from U(0,1)) for N = 2, 6, and 12, and compare histograms to Gaussian curves.
2. **We simulate sums of discrete values** using dice rolls (12-sided) for N = 2, 30, and 50. We again compare results to the normal distribution.
3. **We use a skewed exponential-like distribution** (fX(x) = 0.7e^(-0.7x), x ≥ 0) and simulate sums for N = 5, 40, and 200.
4. **We simulate sums of Bernoulli trials** with p = 0.6 for N = 4, 20, and 100. We compare both the exact probability mass function and the CLT approximation side-by-side in subplots.

For every experiment:
- We generate a large number of trials (e.g., 100,000).
- We compute the sample mean and variance of the sums.
- We compare those with the theoretical values.
- We plot histograms scaled as probability density or mass functions.
- We overlay a matching Gaussian curve with the same mean and variance to show the effect of CLT.
