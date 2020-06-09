# Protocolo Filogenia ALAD

Diego G. Teixeira<br/>
diego.go.tex@gmail.com

Emanuel Duarte <br/>
emmanuel.baduarte@gmail.com

--- 

## Introdução
Esse estudo busca analisar as relações filogenéticas entre a proteína porfobilinogênio sintase (também conhecida como **Delta-aminolevulinic acid dehydratase - ALAD**). Com o objetivo de:

1. Inferir as relações filogenéticas entre as sequências de ALAD de diferentes organismos;
2. Observar as diferenças nos animoácidos dos diferentes grupos;
3. Relacionar os AA consenso com os valores de atração com os co-fatores;  

## Pegando as sequências
O primeiro passo para inferir as relações filogenéticas entre as sequências de proteína da ALAD em diferentes oganismos é baixar as sequências de interesse a partir de um banco de dados. Para isso, utilizei o banco de dados de proteínas do [NCBI](https://www.ncbi.nlm.nih.gov/protein) por "aminolevulinic dehydratase". Quando o NCBI retornar a lista de resultados, eu selecionei a opção **RefSeq** (banco de dados do NCBI com sem sequências repetidas do mesmo organismo, não excluindo diferentes isoformas da mesma proteína). Para baixar as proteínas selecionadas na busca, vá em Send to > File, selecione a opção "GenPept (full)" e clique em "Create File". O formato GenPept é similar ao GenBank.

![NCBI-Protein-RefSeq](https://github.com/diegogotex/filogenia_ALAD/blob/master/imagens/Fig1.png)

## Formatando as sequências
Para transformar o arquivo GenPept em algo que eu possa trabalhar, irei executar algumas linhas de comando em python para formatar em um arquivo tabular contendo algumas informações de interesse. As etapas a seguir foram realzadas em python 3.7 na plafaforma Colab. 

O script vai:
- carregar o arquivo GenPept;
- Puxar as informações de ID, Nome do organismo, TaxonID e Sequêcia proteíca;
- Remover sequências do mesmo organismo (tirar isoformas).


```python

#carregando os arquivos do drive
from google.colab import drive
drive.mount('/content/gdrive')

#instalando o biopython
!pip install biopython

#carregando as bibliotecas
from Bio import SeqIO
import pandas as pd
import numpy as np

#navegando ate a pasta do ALAD
%cd /content/gdrive/My Drive/paper_2/Filogenia_ALAD/seqs

#carregando os arquivi do genbank
gbFile = "sequence_pept.gp"



#loop de teste para saber se deu certo carregar o arquivo completo
count = 0
for seq_record in SeqIO.parse(gbFile, "genbank"):
    #print(seq_record.id)
    count +=1
    #print(count)
    
print("%s sequências" % (count))



#criando um df do pandas
df = pd.DataFrame(data=None, columns=['ID','Organismo','TaxonID','Seq'])

#loop para popular o dataframe
for seq_record in SeqIO.parse(gbFile, "genbank"): 
    count += 1
    #puxando as informações
    source = seq_record.features[0] #puxando as features de cada ACCESSION
    organism = str(source.qualifiers.get('organism')) #pega o valor com informação do organismo
    organism = organism.replace(" ", "_") #removendo os espaços
    organism = str(organism)[2:-2]#removendo a aspa e o colchete
    db_xref = str(source.qualifiers.get('db_xref')) #puxando o taxonID
    db_xref = db_xref.replace("taxon:", "") #removendo os characteres "taxon:"
    db_xref = str(db_xref)[2:-2]#removendo a aspa e o colchete
    #escrevendo no df
    df.loc[count] = [seq_record.id, organism, db_xref, str(seq_record.seq)]



#reordenando o df
df = df.sort_values(by='TaxonID', ascending=True)

#removendo os organismos duplicados
df = df.drop_duplicates(subset='TaxonID', keep='first')

#checando se o arquivo tá certo
df

#escrevendo o DF para um csv
df.to_csv('sequence.csv', index=False)
```

<br/>
No final do processo terá um arquivo .csv o qual eu costumo abrir o libreoffice. O arquivo vai conter a informação de 13678 sequências de diferentes organismos.

**obs:**  nos identificadores XP_001483075.1 e XP_008710822.1 o campo TaxonID tem um erro onde identificador do taxon vem junto a outras informações: "ATCC:6260', '294746" e "HMP:1541', '1220924". Como são só duas sequências eu editei manualment no libreoffice. 

O arquivo .csv eu transformei em .fasta para a etapa de alinhamento.

## Alinhamento
Nessa etapa eu vou alinhas todas as 13678 sequências utilizando o MAFFT e dpeois vou editar o alinhamento com o CIAlign. A edição tem como objetivo remover colunas com um alto volume de GAPs e cotar as extremidades mal alinhadas, deixando o alinhamento mais limpo e fácil de visualizar/trabalhar. 

Os passos a seguir foram realizados no NPAD da UFRN por causa da alta demanda computacional. 

Alinhando as sequências com o mafft:
```bash
#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH --output=mafft_ALAD.out
#SBATCH --error=mafft_ALAD.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64000
#SBATCH --time=10-0:0

module load softwares/mafft/7.397-gnu-7.1

mafft --reorder --thread 32  --auto sequence.fasta > sequence_MAFFT.fasta
```

O alinhamento visualizado no AliView fica dessa forma:

![Alinhamento-mafft-Aliview](https://github.com/diegogotex/filogenia_ALAD/blob/master/imagens/Fig2.png)

Como se pode ver, existem muitas colunas com gaps e blocos com alinhamento pobre. Para melhorar o alinhamento nessa região eu vou utilizar o CIAlign para remover algumas das colunas com pouca quantidade de resíduos, as quais não apresentam nenhuma informação filogenética. 


editando o alinhamento com o trimal:

```bash
#!/bin/bash
#SBATCH --job-name=CIAlign
#SBATCH --output=CIAlign_ALAD.out
#SBATCH --error=CIAlign_ALAD.err
#SBATCH --mem=64000
#SBATCH --time=10-0:0

module load softwares/python/3.6.1-gnu-4.8

~/.local/bin/CIAlign --infile sequence_MAFFT.fasta --outfile_stem sequence_MAFFT_TRIM --remove_insertions --crop_ends

```

![Alinhamento-mafft-CIAlign-Aliview]()

## Inferinfo a filogenia

Para iferir a filogenia dos 