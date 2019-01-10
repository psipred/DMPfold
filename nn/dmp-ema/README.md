## Example usage

`python pytorch_dmpema_pred.py T1010.21c T1010.map`

Outputs model quality class probabilities as follows:

TM < 0.3
0.3 <= TM < 0.4
0.4 <= TM < 0.5
0.5 <= TM < 0.6
0.6 <= TM < 0.7
TM >= 0.7

Finally outputs p(TM >= 0.5) i.e. sum of last 3 bins.
