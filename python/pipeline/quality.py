import pandas as pd, numpy as np
d=pd.read_parquet('/home/claude/work/br/tb_br.parquet')
d['ano']=pd.to_datetime(d.DT_DIAG,errors='coerce').dt.year
d['lab']=(d.BACILOSC_E=='1')|(d.CULTURA_ES=='1')|d.TEST_MOLEC.isin(['1','2'])
d['uf_res']=d.ID_MN_RESI.fillna('').str[:2]
UF={'11':'RO','12':'AC','13':'AM','14':'RR','15':'PA','16':'AP','17':'TO','21':'MA','22':'PI','23':'CE','24':'RN','25':'PB','26':'PE','27':'AL','28':'SE','29':'BA','31':'MG','32':'ES','33':'RJ','35':'SP','41':'PR','42':'SC','43':'RS','50':'MS','51':'MT','52':'GO','53':'DF'}
REG={'1':'Norte','2':'Nordeste','3':'Sudeste','4':'Sul','5':'Centro-Oeste'}
d['uf']=d.uf_res.map(UF); d['regiao']=d.uf_res.str[0].map(REG)
d['se']=pd.to_numeric(d.SITUA_ENCE,errors='coerce')
base=d[d.ano.between(2015,2025)&d.uf.notna()]
cn=base[(base.FORMA=='1')&base.lab&(base.TRATAMENTO=='1')].copy()
print('casos novos pulm lab BR', len(cn), 'com mud diag', (cn.se==6).sum())
print('flow', len(d), len(d[d.ano.between(2015,2025)]), len(base), len(base[base.FORMA=='1']), len(base[(base.FORMA=='1')&base.lab]), len(cn))
def q(g):
    coh=g[g.ano<=2024]
    pre=g[g.ano<=2019].groupby('ano').size().mean(); y20=(g.ano==2020).sum()
    return pd.Series({
     'casos':len(g),
     'queda2020_%':(y20/pre-1)*100,
     'var_2025_vs_pre_%':((g.ano==2025).sum()/pre-1)*100,
     'sem_encerr_%':coh.se.isna().mean()*100 + coh.se.isin([0]).mean()*100,
     'transf_%':(coh.se==5).mean()*100,
     'mud_diag_%':(coh.se==6).mean()*100,
     'raca_ign_%':(~g.CS_RACA.isin(['1','2','3','4','5'])).mean()*100,
     'escol_ign_%':(~g.CS_ESCOL_N.isin([str(i) for i in range(0,9)]+['10'])).mean()*100,
     'hiv_result_%':g.HIV.isin(['1','2']).mean()*100,
     'ppl_ign_%':(~g.POP_LIBER.isin(['1','2'])).mean()*100,
     'ppl_%':(g.POP_LIBER=='1').mean()*100,
     'rua_%':(g.POP_RUA=='1').mean()*100,
     'trm_%':g.TEST_MOLEC.isin(['1','2','3','4']).mean()*100,
     'cultura_%':g.CULTURA_ES.isin(['1','2','3']).mean()*100,
     'cura_%':(coh.se==1).sum()/coh.se.between(1,10).sum()*100,
     'abandono_%':coh.se.isin([2,10]).sum()/coh.se.between(1,10).sum()*100,
     'notif_fora_uf_%':(g.SG_UF_NOT!=g.uf_res).mean()*100,
     'escol_cod_2_4_6_%':g.CS_ESCOL_N.isin(['2','4','6']).mean()*100,
    })
Q=cn.groupby('uf').apply(q)
Q.loc['Brasil']=q(cn)
order=['RO','AC','AM','RR','PA','AP','TO','MA','PI','CE','RN','PB','PE','AL','SE','BA','MG','ES','RJ','SP','PR','SC','RS','MS','MT','GO','DF','Brasil']
Q=Q.loc[order]
pd.set_option('display.width',250,'display.max_columns',30)
print(Q.round(1))
Q.to_csv('/home/claude/work/br/qualidade_uf.csv')
cn.to_parquet('/home/claude/work/br/cn_br.parquet')
print('SITUA codes', d.SITUA_ENCE.value_counts(dropna=False).head(20).to_dict())
