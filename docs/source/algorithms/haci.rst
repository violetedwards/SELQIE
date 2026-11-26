HACI
====

Heat-Bath Adaptive Configuration Interaction (HACI) first screens unique
external determinants using the HCI criterion. A candidate :math:`a` is
retained for ranking if at least one selected determinant :math:`i` satisfies

.. math::

   |H_{ai}c_i| > \epsilon.

The candidates are then ranked by their estimated energy lowering,

.. math::

   I_a =
   \left|\frac{\Delta_a}{2}
   - \sqrt{\frac{\Delta_a^2}{4} + V_a^2}\right|,
   \qquad
   \Delta_a = H_{aa} - E_{\mathrm{var}},
   \qquad
   V_a = \sum_{i \in \mathrm{selected}} H_{ai}c_i.

Candidates are added in descending :math:`I_a` order until the total
importance of the unselected candidates is no greater than :math:`\sigma`.
The existing selected space is always retained.

The screening threshold is set internally to :math:`\epsilon = \sigma`.
For Python, set ``sel_method="haci"`` and provide ``aci_sigma``.
