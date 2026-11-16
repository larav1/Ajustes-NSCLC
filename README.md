## Arquivos

* **`shapes.csv`**: registra as dimensões matriciais (`shape`) e a correspondência geométrica (tamanho, espaçamento, origem e orientação) dos pares de tomografia computadorizada (TC) e máscara tumoral (`SEG`) nas etapas intermediárias de alinhamento (`nii_aligned`) e reamostragem isotrópica (`resampled`).
* **`estatisticas_npy.csv`**: consolida os metadados técnicos dos tensores finais (`ct_multi.npy` e `seg.npy`), incluindo tipos de dados (`float32` e `uint8`), dimensões padronizadas (`128 × 128 × 128`), faixa de intensidades normalizadas (`[0.0, 1.0]`) nos canais de janelamento (tecido mole e pulmonar) e contagem de voxels positivos na máscara.
* **`audit_report.csv`**: apresenta o diagnóstico completo de integridade espacial por paciente, comparando o volume da região de interesse (ROI) antes e após o recorte (*cropping*) e sinalizando casos atípicos (como ROIs vazias, truncamento de lesões nas bordas da janela de recorte e discrepâncias volumétricas).

```
  {
  "total_pacientes_auditados": 421,
  "pacientes_com_ct_e_seg": 421,
  "alinhamento_geometrico_ok": 421,
  "reamostragem_1mm_ok": 421,
  "com_volumes_finais_npy": 421,
  "com_roi_vazia_erro": 52,
  "lesao_totalmente_contida_no_crop": 36,
  "shape_fixo_128_ok": 421,
  "intensidade_normalizada_ok": 421,
  "total_casos_atipicos": 387
}
```
