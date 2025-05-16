# 🔁 CMPE320 Project 5 – Autocorrelation in Random Processes

## 📌 About

In this project, we explore the idea of **autocorrelation**, which is all about seeing how a signal relates to a shifted version of itself over time. We’re working with random signals—specifically white Gaussian noise—and looking at how consistent the patterns are when we shift them left or right.

We use a provided MATLAB script (no coding needed!) to generate a bunch of random signals and calculate their autocorrelation using the `xcorr` function. Then, we look at the results to see if the process is **wide-sense stationary (WSS)** or **ergodic**, and how averaging across time or across multiple random signals affects what we see.

In the second part of the project, we simulate what happens when we **smooth the noise using a sliding window**. This mimics a simple filter and shows how the autocorrelation changes as we increase the size of the window.

---

## 🎯 Goal

The goal is to:
- Understand what the autocorrelation function represents and how it relates to time shifts in random signals.
- See how averaging across multiple random signals gives us a better approximation of the "true" statistical behavior.
- Explore how filtering a random process (using a sliding window average) affects the signal and its autocorrelation.
- Practice interpreting plots and making meaningful observations, even when the math is complex.
- (Optional) Derive the analytical form of the autocorrelation for a filtered process and compare it to what we see in the plots.

---

## ⚙️ How It Works

Here’s how we tackle the project:
1. **We use MATLAB’s `randn`** to generate random Gaussian signals (white noise).
2. **We calculate the autocorrelation** using MATLAB’s `xcorr`, which gives us an estimate of how much the signal “lines up with itself” as we shift it.
3. **We repeat this multiple times**, store the results, and average them to get a clearer picture.
4. **We analyze four main plots** that show how both the original signals and their autocorrelations behave.
5. In the second part, **we apply a sliding window average** to the white noise to simulate filtering. We explore how this smooths the signal and changes the autocorrelation shape.
6. For each filter size (L = 16, 25, 64), we:
   - Observe how the shape of the curve changes.
   - Calculate the variance reduction.
   - Suggest a general rule about filtering and variance.

