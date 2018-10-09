## DADA2 16S Workflow

1. Make a directory and put the de-multiplexed gzipped forward and reverse reads in a subdirectory called raw t
2. Run `python DADA2_prep.py --raw_dir DIR --out_dir OUT --primer_set [amoA nxrB or 16S_515F_926R]`
  * This will trim the adapters from the reads, using whichever primer set specified.
  * You can also pass `--fwd --rev` with the primers you want to use if they are not in the defaults.
3. Run the cells of the DADA2_notebook
  * See the DADA2 tutorial [here](https://benjjneb.github.io/dada2/tutorial.html)
  * You'll have to manually set the `path` var in R (Cell 1:Line 18) to the top-level directory
  * Look at a few of the quality profiles. In the filter step, truncate forward and reverse reads where the quality drops off.
  * In addition to truncating at a given quality score, you can filter based on expected errors. I like to keep this ~1% expected errors, but feel free to change this.
  * Learning the error rates takes forever. Get coffee.
  * The next steps run fast and with no user input. Check `track` at the end to make sure a good fraction of reads are getting through the filtering process. If none of your reads survive a step, check things like trunc-length vs. read length, quality threshold, etc.
4. I use DADA2_biom_fixer.ipynb to separate the ASV reference sequences and make a qiime-compatible biom txt file
5. Use the **biom** utilities to convert the output file to a .biom

### Taxonomy assignment
* Names are cool and sometimes meaningful!
* Two things to consider are reference database and classification method.
* **Reference databases**: In general, I've found recent **Silva** releases > **RDP**
* For methods, **uclust** is the default in QIIME and the creator of uclust thinks his method is superior to the RDP naive bayes [here](https://www.drive5.com/usearch/manual/taxonomy_validation.html)
* DADA2 implements the RDP-Bayes classifier which you're free to use as well. I typically just run `assign_taxonomy.py` in QIIME1 with the default uclust classifier and the most recent Silva database, although this apparently results in less confident (deep taxonomic level) classifications than the RDP-Bayes classifier.
* I may add support for the ASV species assignment capabibility in DADA2 described [here](https://benjjneb.github.io/dada2/assign.html#species-assignment)
