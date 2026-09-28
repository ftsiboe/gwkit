# gwkit 0.0.0.9001

* `approach = "BIC"` is accepted by `estimate_gwss()`, `estimate_gwss_kernels()` and
  `gw_select_bandwidth()` (passed to `GWmodel::bw.gwr()`). Leave-one-out CV can keep
  improving down to the smallest bandwidth when observations sit side by side; BIC
  found an interior optimum for every kernel on elevator basis (arpcPriceBasis,
  2026-09-28) where CV and AICc did not.
* `estimate_gwss_kernels(bw = )` also takes one bandwidth per kernel, as a vector named
  by kernel (e.g. `c(gaussian = 0.6, bisquare = 1.4)`); the per-kernel values are
  attached as the `bandwidth` attribute. A single value is still shared by every kernel.
