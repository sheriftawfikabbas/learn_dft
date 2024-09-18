
*Answers to tutorial attendee questions:*

**Does GPAW support GPU?**

Yes, https://wiki.fysik.dtu.dk/gpaw/platforms/Bull/curie_gpu.html

**What is the default XC in a GPAW DFT calculation?**

It is LDA, as I said in the tutorial.

*Corrections:*

- I mentioned that GPAW supports MD - sorry I meant ASE supports MD.
- In the Pt slab calculation, I quoted the value of H2_energy as -4.48 eV which is wrong, it should be -6.725007250781371 (which was obtained in BasicCalculations_1_H2.py). The value for the adsorption energy makes more sense now: ~ -0.4 eV. However, more k-points are needed to improve it.
