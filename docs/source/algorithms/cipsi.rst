CIPSI
=====

Configuration Interaction using a Perturbative Selection made Iteratively
(CIPSI) generates external determinants through single and double
excitations of the current selected space. For each unique external
determinant :math:`a`, its first-order wavefunction coefficient is

.. math::

   c_a^{(1)} =
   \frac{\sum_{i \in \mathrm{selected}} H_{ai} c_i}
        {E_{\mathrm{var}} - H_{aa}}.

The determinant is added when :math:`|c_a^{(1)}|` exceeds the configured
coefficient threshold. For Python, set ``sel_method="cipsi"`` and provide
``cipsi_coeff_threshold``. A near-zero denominator with nonzero coupling
causes that determinant to be selected as an intruder state.

The iterative perturbative selection scheme was introduced by Huron, Malrieu,
and Rancurel, *J. Chem. Phys.* **58**, 5745 (1973).
