#!/usr/bin/python3.5
import argparse
import json
import sys


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


def _parseCheInfo(genes=[]):

    seqInfo = {
        'neighbors': {
            'chea': [],
            'cheb': [],
            'cher': [],
        },
        'aseqs': [],
        'cheABRac': []
    }

    for gene in genes:
        cheAac = bitk3.getAccessionFromMist22Gene(gene)
        cheABRcounts = {
            'chea': 0,
            'cheb': 0,
            'cher': 0
        }

        for i, cheInfo in enumerate(gene['cheInfo']):
            if cheInfo['che']:
                if cheInfo['che'] == 'cher':
                    accession = gene['neighborsAC'][i]
                    aseq = gene['neighborsAseq'][i]
                    seqInfo['neighbors']['cher'].append({
                        'header': accession,
                        's': aseq,
                        'headCheA': cheAac
                    })
                    cheABRcounts['cher'] += 1
                    seqInfo['aseqs'].append(aseq)
                    seqInfo['cheABRac'].append(accession)
                elif cheInfo['che'] == 'cheb':
                    accession = gene['neighborsAC'][i]
                    aseq = gene['neighborsAseq'][i]
                    seqInfo['neighbors']['cheb'].append({
                        'header': accession,
                        's': aseq,
                        'headCheA': cheAac
                    })
                    cheABRcounts['cheb'] += 1
                    seqInfo['aseqs'].append(aseq)
                    seqInfo['cheABRac'].append(accession)
                elif cheInfo['che'] == 'thisChea':
                    accession = gene['neighborsAC'][i]
                    aseq = gene['neighborsAseq'][i]
                    seqInfo['neighbors']['chea'].append({
                        'header': accession,
                        's': aseq,
                        'headCheA': cheAac
                    })
                    cheABRcounts['chea'] += 1
                    seqInfo['aseqs'].append(aseq)
                    seqInfo['cheABRac'].append(accession)
                elif cheInfo['che'] == 'chea':
                    accession = gene['neighborsAC'][i]
                    aseq = gene['neighborsAseq'][i]
                    seqInfo['neighbors']['chea'].append({
                        'header': accession,
                        's': aseq,
                        'headCheA': cheAac
                    })
                    cheABRcounts['chea'] += 1
                    seqInfo['aseqs'].append(aseq)
                    seqInfo['cheABRac'].append(accession)

        conflict = ''
        for che in cheABRcounts.keys():
            if cheABRcounts[che] == 0:
                seqInfo['neighbors'][che].append({
                    'header': None,
                    's': None,
                    'headCheA': cheAac
                })
            if cheABRcounts[che] > 1:
                conflict += '{}_{}, '.format(cheABRcounts[che], che)

        if conflict != '':
            for che in cheABRcounts.keys():
                # print(che)
                # print(json.dumps(seqInfo['neighbors'][che],indent=2))
                for info in seqInfo['neighbors'][che][-cheABRcounts[che]:]:
                    # print(info)
                    info['header'] += '{}CONFLICT:{}'.format(
                        bitk3.BITKTAGSEP,
                        conflict[:-2]
                    )

    return seqInfo


def _addFastaInfo(seqInfo={}, aseq2seq={}, ac2header={}):
    """ adds FASTA info to SeqInfo """

    seqInfo['FASTA'] = {}

    for che in seqInfo['neighbors'].keys():

        fastaString = ''
        for seq in seqInfo['neighbors'][che]:
            try:
                sequence = aseq2seq[seq['s']]
            except KeyError:
                sequence = 'None'
            if seq['header']:
                ac = seq['header'].split('|CONFLICT')[0]
                fastaString += '>{}\n{}\n'.format(
                    seq['header'].replace(ac, ac2header[ac]) + '::' + seq['headCheA'],
                    sequence
                )
            else:
                fastaString += '>{}\n{}\n'.format(
                    ac2header[seq['headCheA']] + '|NOTFOUND',
                    sequence
                )
        seqInfo['FASTA'][che] = fastaString

    return seqInfo


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
    cheaAC = {}

    client = bitk3.get_mist22_client()
    mist22 = client.mist22

    with open(cheaTagFileName, 'r') as f:
        for tag in f:
            tag = tag.replace('\n', '')
            accession = bitk3.bitk3tagToAccession(tag)
            cheaAC[accession] = tag

    genes, badAC = bitk3.accessionToGeneInfo(cheaAC)
    for gene in genes:
        gene = _getNeighbors(mist22, gene, geneNeighborhoodWindow)
        gene = _getSigTransInfoOfNeighbors(mist22, gene)

    seqInfo = _parseCheInfo(genes)
    aseq2seq = bitk3.getSeqFromAseq(seqInfo['aseqs'])

    cheBRgenes, badAC = bitk3.accessionToGeneInfo(seqInfo['cheABRac'])
    cheBRgenes = bitk3.addBitk3tagTomist22GeneInfo(cheBRgenes)

    ac2bitk3tag = {None: None}
    for i, ac in enumerate(seqInfo['cheABRac']):
        ac2bitk3tag[ac] = [
            g['bitk3tag'] for g in cheBRgenes if g['p']['ac'] == ac
        ][0]

    seqInfo = _addFastaInfo(seqInfo, aseq2seq, ac2bitk3tag)

    for che in seqInfo['neighbors'].keys():    
        filename = che + '.fa'
        with open(filename, 'w') as f:
            f.write(seqInfo['FASTA'][che])

    client.close()

    return seqInfo


if __name__ == "__main__":
    import bitk3
    parser = argparse.ArgumentParser(
        prog='buildCheABRmakeFasta',
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
