## Etapas do pré-processamento

| **Etapa** | **Descrição** | **Resultado esperado** | **Status (421 pacientes)** |
|:---|:---|:---|:---|
| **1. Dados originais** | Utilização da coorte NSCLC-Radiomics (LUNG1), com TC e segmentação tumoral | Cada paciente deve ter os arquivos correspondentes. | **OK** |
| **2. Conversão de formato e leitura dos arquivos** | Leitura dos NIfTI de TC e máscara (`nii_aligned`) | Arquivos legíveis e associação correta ao paciente | **OK** |
| **3. Alinhamento espacial** | TC e máscara são colocadas no mesmo referencial espacial. | A segmentação deve indicar a localização correta do tumor na TC. | **OK** |
| **4. Verificação geométrica** | Comparação de tamanho, origem, espaçamento e orientação. | TC e SEG alinhadas devem apresentar geometria compatível. | **OK** (Geometria idêntica em 421/421) |
| **5. Reamostragem da TC** | Conversão da resolução para $1\text{ mm}^3$ isotrópico (`resampled`). | Preservar a anatomia e a informação de intensidade. | **OK** (Espaçamento `1.0 x 1.0 x 1.0 mm`) |
| **6. Reamostragem da SEG** | Conversão da máscara para a mesma resolução (`resampled`). | Preservar a região tumoral utilizando interpolação *nearest neighbor*. | **OK** (421/421 máscaras íntegras) |
|**7. Clipping e Janelamento** | Global Clipping ($-1200$ a $600$ HU) e Multiple Windowing (janela pulmonar e janela de tecido mole) | Realçar estruturas relevantes sem distorcer os limites físicos das Unidades Hounsfield | OK |
| **8. Normalização**  | Ajuste dos valores de intensidade dos canais. | Valores dentro da faixa prevista, neste caso $[0, 1]$. | OK |
| 9. Crop / Pad ---  | Recorte ou preenchimento dos volumes para dimensão fixa de $128^3$. | Preservar integralmente a região tumoral e aplicar operações espacialmente correspondentes na TC e SEG. | **ERRO** (52 ROIs vazias; 333 tumores com truncamento nas bordas) |
| 10. Construção multicanal ---  | Organização dos canais de TC (janela pulmonar e janela de tecido mole/mediastino). | Canais devem corresponder exatamente ao mesmo volume espacial. | Refazer |
| 10. Preparação da SEG ---  | Organização da máscara tumoral final em $128^3$. | A máscara deve continuar contendo os voxels positivos originais, sem amputação ou perda da lesão. | **COMPROMETIDA PELA ETAPA 7** (52 casos com `0` voxels; 333 parciais) |
| 11. Salvamento --- | Geração de `ct_multi.npy` e `seg.npy` (`processed_final`). | Arquivos corretos, shapes padronizados e associação por paciente. | Refazer |


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

1. com_roi_vazia_erro (52): Em 52 pacientes, a máscara final (seg.npy) terminou com zero voxels de tumor, ou seja, o recorte caiu fora da lesão ou a máscara sumiu.
2. lesao_totalmente_contida_no_crop (36): Apenas 36 pacientes tiveram o tumor inteiro preservado dentro do cubo de $128^3$ sem tocar ou ser cortado pelas bordas.
3. total_casos_atipicos (387): Soma de todos os casos que apresentaram algum desvio nos critérios espaciais ou volumétricos (os 52 vazios + os 333 cortados + 2 casos com alguma outra anomalia de volume extremo)

OBS: apenas 8,6% da coorte (36/421) contém o tumor íntegro sem cortes ou deslocamento. Por isso o modelo de unimodal imagem teve o pior desempenho.

## O que deu errado no código original (`image_preprocessing.ipynb`)

Na etapa de recorte e padronização dimensional (`get_bounding_box` + `crop_to_bbox` + `pad_or_crop_to_shape`), ocorreram dois problemas combinados:

1. Corte central rígido de 128 voxels (`if size >= target`): Quando a *bounding box* do tumor (acrescida da margem de 5 voxels) ultrapassava `128 mm` em qualquer eixo, a função `pad_or_crop_to_shape` aplicava um recorte central fixo (`start = (size - 128) // 2`), removendo as extremidades da lesão (causa dos 333 casos truncados).

```
if size >= target:
    start = (size - target) // 2
    crop_slices.append(slice(start, start + target))
else:
    start = (target - size) // 2
    pad_slices.append(slice(start, start + size))
```

2. Deslocamento do centro em máscaras não convexas ou com múltiplos componentes: Quando a máscara possuía mais de um componente desconexo ou formato irregular, o centro geométrico da *bounding box* caía em uma região vazia entre os componentes, fazendo com que o corte central de `128³` capturasse apenas tecido saudável/fundo e zerasse a máscara (causa das 52 ROIs vazias).

```
def get_bounding_box(mask, margin=5):
    coords = np.array(np.where(mask > 0))
    if coords.size == 0:
        return None

    min_coords = coords.min(axis=1) - margin
    max_coords = coords.max(axis=1) + margin
    ...
```

Bounding Box (Caixa Delimitadora): Envolve completamente todos os voxels positivos da máscara de segmentação tumoral, delimitado pelas coordenadas mínimas e máximas da lesão.

Crop (Recorte): É a operação de extrair apenas a sub-região de interesse (ROI) do volume tridimensional da tomografia e da máscara (definida pela bounding box ou por uma janela centrada no tumor), descartando o restante do fundo anatômico irrelevante (como bordas externas do tórax ou regiões distantes da lesão).

Padding (Preenchimento): É a adição controlada de voxels constantes nas bordas de um volume recortado quando suas dimensões são menores do que o tamanho fixo exigido pela rede neural (por exemplo, $128 \times 128 \times 128$), utilizando valores neutros de fundo (como 0 para máscaras binárias/intensidades normalizadas ou -1024 HU para representar o ar na tomografia).

## Sugestão de ajuste 

- [ ] **Evitar o corte nas bordas:** Quando a *bounding box* do tumor (+ margem) for maior que `128³`, redimensionar (*resize*) o recorte para `(128, 128, 128)` em vez de cortar as extremidades, usando interpolação linear na TC e *nearest neighbor* na máscara.
- [ ] **Evitar ROI vazia:** Garantir que o centro do recorte seja calculado a partir dos voxels positivos da máscara (`mask > 0`), conferindo no final que `seg_final.sum() > 0` para todos os pacientes. Inserir uma verificação automática para identificar pacientes cuja máscara final fique vazia `(seg_final.sum() == 0)`.
- [ ] **Validar visualmente:** Plotar 3 vistas (axial, coronal e sagital) de alguns casos que antes davam erro (ex.: `LUNG1-007`) para confirmar que o tumor está centralizado e inteiro.
