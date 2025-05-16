# 📡 CMPE320 Project 4 – MAP and ML Detection of Binary Signals

## 📌 About

In this project, we explore how to detect binary data that’s been transmitted through a noisy channel using techniques from probability and communication theory. We simulate a digital communication system where bits (0s and 1s) are sent by encoding them as signal levels (+A or -A), but noise (random fluctuations) gets added to these signals along the way. Our job is to build smart detectors that try to recover the original bits as accurately as possible from the noisy data.

We design and test two types of detectors:  
- **MAP (Maximum A Posteriori)**: Uses prior probabilities to help decide what bit was most likely sent.  
- **ML (Maximum Likelihood)**: Assumes both bits are equally likely and just chooses the one that best fits the received signal.

We vary the noise level, simulate the transmission and recovery of bits, and measure how often each method gets it right. This gives us a real, hands-on understanding of how bit error rates change depending on the signal-to-noise ratio.

---

## 🎯 Goal

The goal is to:
- Derive and implement the MAP and ML detection rules based on conditional probabilities.
- Visualize how noise affects the shape of signal distributions and where the detectors set their decision thresholds.
- Simulate bit transmission and decoding, and calculate how often errors happen for different noise levels.
- Compare the theoretical and simulated probability of bit errors (BER) across a range of signal-to-noise ratios (SNR).
- Observe when and why MAP performs better than ML, especially when the probability of sending a 0 or 1 isn’t equal.

---

## ⚙️ How It Works

Here’s what we do:
1. **We model the transmission system** where 0 is sent as +A and 1 is sent as -A, with A = 3V. Noise is added using a Gaussian distribution with different variances.
2. **We derive the MAP detection rule** using probability theory, including Bayes’ Rule. We solve for the decision threshold (tau_MAP) that tells us when to guess 0 or 1 based on the received voltage.
3. **We simulate the system at various prior probabilities (p0)** and plot how the MAP threshold changes depending on how likely we are to send a 0.
4. **We test the MAP detector** by running simulations with thousands of trials at a specific SNR (e.g., 10 dB), plot the histogram of received signals, and show the MAP threshold on the plot.
5. **We derive and simulate the ML detector**, which is just a special case of MAP when the bits are equally likely. We derive a clean equation for bit error rate using the Q function and plot both simulated and analytical results across a wide SNR range.
6. **We repeat the process for MAP**, using a non-uniform prior (like p0 = 0.25), and compare its bit error rate against ML’s.
7. **We analyze** how MAP and ML perform differently and plot the ratio of their bit error probabilities to show when MAP gives a real advantage.
