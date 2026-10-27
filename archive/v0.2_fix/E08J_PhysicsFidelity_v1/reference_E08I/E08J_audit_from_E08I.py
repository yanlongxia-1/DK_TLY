"""Retrospective audit of existing E08I CSV, no new data and no solver run.
Method was designed with already examined Z9-1, NOT held-out evidence.
"""
import math
from pathlib import Path
import pandas as pd
BASE=Path(__file__).resolve().parent
src=pd.read_csv(BASE/'E08I_band_detail.csv')
rows=[]; summary=[]
for key,g in src.groupby(['fileName','K','alpha','reconstruction'],sort=True):
    g=g.copy()
    mx=g.groupby('channel')['rawFitAmplitude'].transform('max')
    strong=g.rawFitAmplitude>=.10*mx
    retention=g.retentionPct/100
    theta=g.phaseChangeDeg.map(math.radians)
    e=g.rawFitAmplitude*(((retention*theta.map(math.cos)-1)**2+(retention*theta.map(math.sin))**2)**0.5)
    frac=e/mx
    strong_good=(g.retentionPct.sub(100).abs()<=1) & (g.phaseChangeDeg<=.2)
    good=(strong & strong_good) | (~strong & (frac<=.01))
    paired={f: bool(strong.loc[z.index].all()) for f,z in g.groupby('frequencyHz')}
    phase_good=all(not paired[f] or z.relativePhaseChangeDeg.max()<=.2
                   for f,z in g.groupby('frequencyHz'))
    g['isStrong']=strong
    g['complexErrorOverReference']=frac
    g['revisedCellPass']=good
    rows.append(g)
    summary.append(dict(fileName=key[0],K=key[1],alpha=key[2],reconstruction=key[3],
       oldAllPass=bool(g.allPhysicsPass.all()),
       revisedAllPass=bool(good.all() and phase_good),
       nStrong=int(strong.sum()),nWeak=int((~strong).sum()),
       maxStrongDeviationPct=float(g.loc[strong,'retentionPct'].sub(100).abs().max()) if strong.any() else float('nan'),
       maxWeakComplexFraction=float(frac.loc[~strong].max()) if (~strong).any() else 0.0,
       nOldHighAmpViolations=int(g.highAmpViolation.sum())))
pd.concat(rows).to_csv(BASE/'E08J_E08I_band_retrospective.csv',index=False)
pd.DataFrame(summary).to_csv(BASE/'E08J_E08I_case_retrospective.csv',index=False)
print(pd.DataFrame(summary)[['fileName','K','alpha','reconstruction','oldAllPass','revisedAllPass']].to_string(index=False))
