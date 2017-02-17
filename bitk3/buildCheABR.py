#!/usr/bin/python3.5
import argparse
import json


def _getSigTransInfoOfNeighbors(mist22Client, gene={}):
    # print(gene['neighborsId'])
    sigTransInfo = mist22Client.signal_genes5.find(
        {
            'cid': gene['cid'],
            '_id': {
                '$in': gene['neighborsId']
            }
        }
    )

    cheInfo = [{'che': False} for i in range(len(gene['neighborsId']))]

    for info in sigTransInfo:
        if info['r'][0] == 'chemotaxis':
            if info['r'][1] == 'cheb':
                cheInfo[gene['neighborsId'].index(info['_id'])]['che'] = 'cheb'
            elif info['r'][1] == 'cher':
                cheInfo[gene['neighborsId'].index(info['_id'])]['che'] = 'cher'

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

    ids = []
    aseqs = []
    for neighbor in neighbors:
        # print(json.dumps(neighbor, indent=2))
        _id = bitk3.getMistIDFromMist22Gene(neighbor)
        ids.append(_id)
        aseq = bitk3.getAseqFromMist22Gene(neighbor)
        aseqs.append(aseq)
    gene['neighborsId'] = ids
    gene['neighborsAseq'] = aseqs

    return gene


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
    geneInfo = {}

    client = bitk3.get_mist22_client()
    mist22 = client.mist22

    with open(cheaTagFileName, 'r') as f:
        for tag in f:
            accession = bitk3.bitk3tagToAccession(tag.replace('\n', ''))
            geneInfo[accession] = {'tag': tag}

    tags = [i for i in geneInfo.keys()]
    genes = bitk3.accessionToGeneInfo(tags)

    for gene in genes:
        # print(json.dumps(gene, indent=2))
        aseq = bitk3.getAseqFromMist22Gene(gene)
        cheAseqs.append(aseq)

        accession = bitk3.getAccessionFromMist22Gene(gene)
        geneInfo[accession]['aseqs'] = aseq

        gene = _getNeighbors(mist22, gene, geneNeighborhoodWindow)
        gene = _getSigTransInfoOfNeighbors(mist22, gene)

        # print(json.dumps(gene, indent=2))
        for i, cheInfo in enumerate(gene['cheInfo']):
            if cheInfo['che']:
                aseq = bitk3.getAseqFromMist22Gene(gene)
                cheAseqs.append(aseq)
                print('{} - {} - {}'
                      .format(
                          cheInfo['che'],
                          gene['neighborsAseq'][i],
                          gene['neighborsId'][i]
                      ))

    print('Closing client')
    client.close()

    return 0


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
