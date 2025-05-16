# 📶 CMPE320 Project 2 – Functions of a Random Variable

## 📌 About

In this project, we explore what happens when we apply different functions to a random signal that's been corrupted by noise. We start with a signal that randomly flips between +1.5V and -1.5V. Then we add Gaussian noise to it—this models what would happen in real communication systems where signals get distorted on their way to the receiver. That final signal is what we call R, and it's the input to everything we do in the rest of the project.

We look at three methods of detecting or processing this noisy signal:
1. A **diode detector** that only keeps the positive part of the signal,
2. An **absolute value detector** that measures the amplitude regardless of sign,
3. A **square law detector** that focuses on signal power by squaring the input.

For each method, we treat the output as a new random variable and find its probability distribution. We simulate everything using MATLAB, and we also do some math by hand to get the analytical PDFs. Then we compare what we got from simulation to what we got from theory. In the end, we reflect on how these transformations change the signal, and we explore whether Jensen’s Inequality holds when we compare E[g(R)] to g(E[R]).

---

## 🎯 Goal

The goal is to:
- Understand how applying a function to a random variable changes its distribution.
- Explore three detection methods: a perfect diode detector, an absolute value detector, and a square law detector.
- Derive analytical expressions for the output PDFs using the CDF method.
- Simulate each method in MATLAB and compare the analytical results to histogram-based PDFs.
- Observe the relationship between the expected value of the output and the function applied to the expected value of the input.
- See if Jensen’s Inequality holds across all three cases.

---

## ⚙️ How It Works

Here's how we approach the project:
1. **We simulate the received signal (R)** by combining a random signal (+A or -A) with Gaussian noise.
2. **We visualize the distribution of R** using both scatterplots and histograms, then overlay the theoretical PDF to confirm our model.
3. **We apply Method 1 (Diode Detector):** Output is 2kR if R > 0, otherwise 0.
4. **We apply Method 2 (Absolute Value Detector):** Output is 4|R|.
5. **We apply Method 3 (Square Law Detector):** Output is 3R².
6. For each method:
   - We derive the analytical PDF using the CDF method.
   - We simulate the output values in MATLAB.
   - We plot the histogram of simulated values and overlay the analytical PDF.
   - We calculate the expected value of the output, E[S], and compare it with g(E[R]).
7. **We wrap up** by comparing results across all three methods and observing whether Jensen’s Inequality appears to hold in each case.
