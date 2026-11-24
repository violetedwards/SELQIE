Selection Utilities
===================

This module provides shared utilities for the selected-CI algorithms.

``any_selected_det``
--------------------

Returns true when a determinant is present in a determinant list.

``remove_selected_dets``
------------------------

Removes determinants already present in a selected space from a candidate list.

``sort_dets_by_score``
----------------------

Sorts determinant/score pairs from highest to lowest score using an in-place
heapsort (:math:`O(n \log n)`). The caller chooses whether scores are compared
by their absolute values.
