# CMEPDAexam
Repository on the project I made for the exam of Computing Methods for Experimental Physics and Data Analysis.

The task is to reproduce the measurement of the differential production cross section for $Z^{0}$ in p-p collision at 13 TeV made by CMS with open data. 

For the article I refer to: https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-17-010/index.html.

We select muon pairs applying the cuts described in the article:

- n = 2
- $p_{t_{i}}$ > 25 GeV $i = 1, 2$
- $|\eta_{i}|$ < 2.4 $i = 1, 2$
- $\Delta R_{i}$ < 0.3 $i = 1, 2$
- $q_{1} + q_{2}$ = 0
  

After finding the resonance of the $Z^{0}$, we calculate distributions for the following quantities:

- transverse momentum : $P_{t}$
-  rapidity : $y$
- optimized angle : $\phi_{\eta} tan\left(\frac{\pi - \Delta \phi}{2}\right)sin\left(\theta_{\eta}\right)$

With $cos\left(\theta_{\eta}\right) = tanh\left(\frac{\Delta \eta}{2}\right)$.

The differential cross section is obtain by scaling the histograms with the integrated luminosity. To reduce inefficiencies due to reconstruction algorithm we performed an unfolding procedure, estimating the response matrix with data generated from a Montecarlo and using the RooUnfold toolbox.

