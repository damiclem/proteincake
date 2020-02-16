import os
import argparse

def test_architecture(module, pfam_id, in_path, background, out, out_wordcloud, p_value = 0.05, depth = 4):
    TARGET =  '{}/{}.csv'.format(in_path, pfam_id)
    OUT = '{}/result_{}.csv'.format(out, pfam_id)
    OUT_WORDCLOUD = '{}/worldcloud_{}.png'.format(out_wordcloud, pfam_id)

    os.system("python {} --out_path {} --out_wordcloud {} --target_path {} --background_path {} --p_value {} --depth {}".
              format(module, OUT, OUT_WORDCLOUD, TARGET, background, p_value, depth))

if __name__ == '__main__':

    # 1. Define arguments
    parser = argparse.ArgumentParser()
    ### Enrichment module to call
    parser.add_argument('--module',               type=str,   default='modules/do_modules/enrichment_do.py')
    ### In directory
    parser.add_argument('--in_path',              type=str,   default='data/architecture/do_architectures')
    parser.add_argument('--background',           type=str,   default='data/architecture/do_architectures/architecture_background.csv')
    ### Out directory of the results and of the WordCloud
    parser.add_argument('--out_path',             type=str,   default='results/do_enrichment/architectures')
    parser.add_argument('--out_wordcloud',        type=str,   default='results/do_enrichment/architectures')
    ### Parameters of the filter
    parser.add_argument('--p_value',              type=float, default=0.05)
    parser.add_argument('--depth',                type=int,   default=4)
    ### Names of the columns of the input dataframe containing the GO_id and the description


    # 2. Define dictionary of args
    args = parser.parse_args()
    datasets = os.listdir(args.in_path)
    for n, dataset in enumerate(datasets):
        if dataset[0] == 'P':
            print('STARTING ITERATION {}/{} on architecture {}'.format(n, len(datasets)-1, dataset))
            test_architecture(module=args.module, pfam_id=dataset[:-4],
                              in_path=args.in_path, background=args.background,
                              out=args.out_path, out_wordcloud=args.out_wordcloud,
                              p_value = args.p_value, depth=args.depth)
