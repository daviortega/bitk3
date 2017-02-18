#!/usr/bin/python3.5
import argparse
import SeqDepot

import json


def _getSigTransInfoOfNeighbors(mist22Client, gene={}):
    # print(gene['neighborsId'])
    sigTransInfo = mist22Client.signal_genes6.find(
        {
            'cid': gene['cid'],
            '_id': {
                '$in': gene['neighborsId']
            }
        }
    )

    cheInfo = [{'che': False} for _id in gene['neighborsId']]

    for info in sigTransInfo:
        if info['r'][0] == 'chemotaxis':
            if info['r'][1] == 'cheb':
                cheInfo[gene['neighborsId'].index(info['_id'])]['che'] = 'cheb'
            elif info['r'][1] == 'cher':
                cheInfo[gene['neighborsId'].index(info['_id'])]['che'] = 'cher'
            elif info['_id'] == gene['_id']:
                cheInfo[gene['neighborsId'].index(info['_id'])]['che'] = 'thisChea'
            elif info['r'][1] == 'chea':
                cheInfo[gene['neighborsId'].index(info['_id'])]['che'] = 'chea'

    gene['cheInfo'] = cheInfo

    return gene


def _getNeighbors(mist22Client, gene={}, geneNeighborhoodWindow=5):
    startSearch = gene['pos'] - geneNeighborhoodWindow
    limitSearch = 2 * geneNeighborhoodWindow + 1
    componentID = gene['cid']
    neighbors = mist22Client.genes.find(
        {
            'cid': componentID,
            'pos': {
                '$gte': startSearch
            }
        }
    ).limit(limitSearch)

    accessions = []
    aseqs = []
    ids = []
    for neighbor in neighbors:
        # print(json.dumps(neighbor, indent=2))
        accession = bitk3.getAccessionFromMist22Gene(neighbor)
        accessions.append(accession)
        aseq = bitk3.getAseqFromMist22Gene(neighbor)
        aseqs.append(aseq)
        _id = bitk3.getMistIDFromMist22Gene(neighbor)
        ids.append(_id)
    gene['neighborsAC'] = accessions
    gene['neighborsAseq'] = aseqs
    gene['neighborsId'] = ids

    return gene


def _parseGeneInfo(genes={}):

    seqInfo = {
        'neighbors': {
            'chea': [],
            'cheb': [],
            'cher': [],
        },
        'aseqs': [],

    }

    cheBRac = []

    for gene in genes:
        aseq = bitk3.getAseqFromMist22Gene(gene)
        cheAseqs.append(aseq)
        cheAac = bitk3.getAccessionFromMist22Gene(gene)
        seqInfo['chea'].append({
            'header': cheAac,
            's': aseq}
        )

        # print('getting info')
        gene = _getNeighbors(mist22, gene, geneNeighborhoodWindow)
        gene = _getSigTransInfoOfNeighbors(mist22, gene)
    #        print(gene['cheInfo'])

        for i, cheInfo in enumerate(gene['cheInfo']):
            if cheInfo['che']:
                aseq = bitk3.getAseqFromMist22Gene(gene)
                cheAseqs.append(aseq)
                accession = gene['neighborsAC'][i]
                cheBRac.append(accession)
                if cheInfo['che'] == 'cher':
                    seqInfo['cher'].append({
                        'header': accession,
                        's': aseq
                    })
                if cheInfo['che'] == 'cheb':
                    seqInfo['cheb'].append({
                        'header': accession,
                        's': aseq
                    })

        if len(seqInfo['chea']) != len(seqInfo['cheb']):
            seqInfo['cheb'].append({
                'header': cheAac,
                's': None
            })
        if len(seqInfo['chea']) != len(seqInfo['cher']):
            seqInfo['cher'].append({
                'header': cheAac,
                's': None
            })

    return SeqInfo


def main(cheaTagFileName='', geneNeighborhoodWindow=5):
    """
    It will build a concatenated alignment of CheABR

    Argument Keywords:
    cheaTagList -- File containing list of bitk3 tags, 1 per line
    geneNeighborhoodWindow -- Integer of the gene neighborhood window to \
    search for CheB or CheR

    Output:
    File with concatenated alignment of ABR
    JSON file with information

    Return:
    0 is it does not fail

    """

    cheAseqs = []
    ac2bitk3tag = {}

    client = bitk3.get_mist22_client()
    mist22 = client.mist22

    with open(cheaTagFileName, 'r') as f:
        for tag in f:
            tag = tag.replace('\n', '')
            accession = bitk3.bitk3tagToAccession(tag)
            ac2bitk3tag[accession] = tag

    cheAac = [i for i in ac2bitk3tag.keys()]
    genes, badAC = bitk3.accessionToGeneInfo(cheAac)

    seqInfo = {
        'chea': [],
        'cheb': [],
        'cher': [],
    }

    cheBRac = []

    for gene in genes:
        aseq = bitk3.getAseqFromMist22Gene(gene)
        cheAseqs.append(aseq)
        cheAac = bitk3.getAccessionFromMist22Gene(gene)
        seqInfo['chea'].append({
            'header': cheAac,
            's': aseq}
        )

        # print('getting info')
        gene = _getNeighbors(mist22, gene, geneNeighborhoodWindow)
        gene = _getSigTransInfoOfNeighbors(mist22, gene)
#        print(gene['cheInfo'])

        for i, cheInfo in enumerate(gene['cheInfo']):
            if cheInfo['che']:
                aseq = bitk3.getAseqFromMist22Gene(gene)
                cheAseqs.append(aseq)
                accession = gene['neighborsAC'][i]
                cheBRac.append(accession)
                if cheInfo['che'] == 'cher':
                    seqInfo['cher'].append({
                        'header': accession,
                        's': aseq
                    })
                if cheInfo['che'] == 'cheb':
                    seqInfo['cheb'].append({
                        'header': accession,
                        's': aseq
                    })

        if len(seqInfo['chea']) != len(seqInfo['cheb']):
            seqInfo['cheb'].append({
                'header': cheAac,
                's': None
            })
        if len(seqInfo['chea']) != len(seqInfo['cher']):
            seqInfo['cher'].append({
                'header': cheAac,
                's': None
            })

        print('{}-{}-{}'.format(len(seqInfo['chea']), len(seqInfo['cheb']), len(seqInfo['cher'])))

    aseq2seq = {}
    sd = SeqDepot.new()
    seqs = sd.find(cheAseqs, {'fields': 's'})
    if seqs:
        for seq in seqs:
            aseq2seq[seq['data']['id']] = seq['data']['s']

    cheBRgenes, badAC = bitk3.accessionToGeneInfo(cheBRac)
    cheBRgenes = bitk3.addBitk3tagTomist22GeneInfo(cheBRgenes)

    for i, ac in enumerate(cheBRac):
        ac2bitk3tag[ac] = [g['bitk3tag'] for g in cheBRgenes if g['p']['ac'] == ac][0]

    print(json.dumps(seqInfo, indent=2))

    for che in seqInfo.keys():
        filename = che + '.fa'
        fastaString = ''
        for seq in seqInfo[che]:
            try:
                sequence = aseq2seq[seq['s']]
            except KeyError:
                sequence = 'None'
            fastaString += '>{}\n{}\n'.format(
                ac2bitk3tag[seq['header']],
                sequence
            )
        with open(filename, 'w') as f:
            f.write(fastaString)

    client.close()

    return SeqInfo


if __name__ == "__main__":
    import bitk3
    parser = argparse.ArgumentParser(
        prog='buildCheABR',
        usage='%(prog)s input_file position',
        description='Makes concatenated alignment of CheABR \
        given a list of CheA tags'
    )
    parser.add_argument(
        'cheaTagFileName',
        metavar='myCheATags.txt',
        type=str,
        help='1 cheA bitk3-tag per line'
    )
    parser.add_argument(
        'geneNeighborhoodWindow',
        metavar='5',
        type=int,
        nargs='?',
        default=5,
        help='Integer of the gene neighborhood window \
        to search for CheB or CheR'
    )

    args = parser.parse_args()

    cheaTagFileName = args.cheaTagFileName
    geneNeighborhoodWindow = args.geneNeighborhoodWindow

    main(cheaTagFileName, geneNeighborhoodWindow)
else:
    from bitk3 import bitk3
