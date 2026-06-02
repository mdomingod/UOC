# Training Run Report — one_run_baseline

_Generated: 2026-05-26 21:00:38_

## Training

| Item | Value |
|---|---|
| Epochs run        | 50 |
| Training time     | 1337.7 s |
| Best val MSE      | 20,135.66 |
| Train samples     | 146,518 |
| Val samples       | 37,370 |

## Site-level metrics

How well does the model predict individual UMI counts?

| Metric | Value |
|---|---|
| n_sites     | 37,370 |
| MSE         | 20,138.84 |
| RMSE        | 141.91 |
| MAE         | 56.95 |
| Pearson r   | 0.6334 |
| Spearman r  | 0.5998 |

## Primer-level metrics

Can the model rank primers by off-target risk?

| Metric | Value |
|---|---|
| n_primers              | 966 |
| n_dropped (zero total) | 0 |
| Mean bias              | +0.0131 |
| Mean abs bias          | 0.0400 |
| RMSE on rates          | 0.0899 |
| Pearson r              | 0.1872 |
| Spearman r             | 0.1708 |
