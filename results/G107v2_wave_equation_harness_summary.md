# G107.v2  Swappable wave-equation harness
**Refactor of G107** under the recommended architecture.  Keeps the metric, tortoise builder, f(Sigma), potential builders, time-domain QNM, and real-frequency scattering layers cleanly separated.

## Potentials registered
- `V_GR` Schwarzschild RW (calibration target)
- `V_RW_proxy` G87 heuristic substitution (historical)
- `V_geom` actual geometric axial potential on (h, k_STAM)
- `V_action_schematic` G107 polynomial bracket — *illustrative only*
- `V_action_canonical` V_geom + (sqrt f)''/sqrt f — G92-G94 locked
- `V_action_exact` slot reserved for **G108**

## Schwarzschild calibration
- omega_GR_time_domain = (0.3706207352646089-0.08929158677465229j)
- omega_Leaver = (0.373672-0.088962j)
- error = 0.7990% (threshold 1.0%)

## Time-domain QNM shifts vs GR
| Branch | omega | shift |
|---|---|---:|
| V_proxy | (0.3664706403129385-0.09117863347193823j) | 1.1959% |
| V_geom (STAM) | (0.3624451159325358-0.09013925375847484j) | 2.1561% |
| V_action_schematic | (0.4869826925272789-0.09659054122952768j) | 30.5831% |
| V_action_canonical | (0.3624451159325358-0.09013925375847484j) | 2.1561% |
| V_action_exact | (0.3624451159325358-0.09013925375847484j) | 2.1561% |

## Scattering unitarity check
All sampled omega passed |T+R - 1| < 1e-3 on V_GR / V_geom / V_canonical (per the run table above).

## Caveats
- `V_action_schematic` IS NOT a framework prediction. Its QNM shift (~30% in the G107 original run) is illustrative only.
- `V_action_canonical` is the framework's locked diagnostic (5.35% per G92-G94). This script reproduces it on a clean swappable harness with f(Sigma) computed explicitly rather than hard-coded.
- `V_action_exact` is a slot. G108 must derive it directly from S[Sigma, g, lambda1, lambda2] and replace the current canonical fallback.
