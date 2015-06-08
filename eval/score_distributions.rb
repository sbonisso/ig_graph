require_relative 'tools/run_iggraph'
require 'tempfile'
require 'bio'
require 'csv'
require 'fileutils'
#
# compute score distributions for given file
# output as CSV with cols as: [reads, rand_dist, shuff_dist]
# where rand_dist is a sequence of same length from same nuc distribution
# and shuff_dist is shuffled sequences
#
def compute_score_dists(read_f, organism, outfile="output_vdj_scores.csv")
  #
  # create rand_dist output
  rand_f = Tempfile.new("rand")
  i = 0
  Bio::FlatFile.open(Bio::FastaFormat, read_f) do |ff|
    ff.each do |entry|
      seq = entry.seq
      rseq = resample_rand_read(seq)
      rand_f.puts ">#{i}\n#{rseq}"
      i += 1
    end
  end
  rand_f.close  
  #
  # create rand_dist output
  shuff_f = Tempfile.new("rand")
  i = 0
  Bio::FlatFile.open(Bio::FastaFormat, read_f) do |ff|
    ff.each do |entry|
      seq = entry.seq
      shuff_seq = seq.split("").shuffle.join("")
      shuff_f.puts ">#{i}\n#{shuff_seq}"
      i+=1
    end
  end
  shuff_f.close
  #
  # output of reads
  out_vdj = run_iggraph(read_f, organism)
  out_rand = run_iggraph(rand_f.path, organism)
  out_shuff = run_iggraph(shuff_f.path, organism)
  #
  # output distributions of V/D/J scores
  CSV.open(outfile, "w") do |csv|
    csv << ["index", "V", "D", "J", "type"]
    lbl = ["read", "rand", "shuff"]
    [out_vdj, out_rand, out_shuff].each_with_index do |out_csv,i|
      index = 0
      #[out_vdj].each do |out_csv|
      CSV.foreach(out_csv.path, {:col_sep => "\t"}) do |row|
        csv << [index, row[-3].to_f, row[-2].to_f, row[-1].to_f, lbl[i]]
        index += 1
      end
    end
  end
end
#
#
#
def resample_rand_read(seq)
  len = seq.size
  
  h = get_prob_h(seq)
  #puts h.to_s
  r = Random.new
  rseq = (0..len-1).map do |i|
    rv = r.rand
    ch = "A"
    s = 0
    [:A, :C, :G, :T].each do |c|
      r1 = (s..s+h[c])
      if r1 === rv then
        ch = c.to_s
        break
      else
        s += h[c]
      end
    end
    ch
  end.join("") 
  rseq
end
#
def get_prob_h(seq)
  h = {:A => 0, :C => 0, :G => 0, :T => 0}
  seq.each_char do |ch|    
    next if h[ch.to_sym].nil?
    h[ch.to_sym] += 1
  end
  total = h.values.inject(0){|s,v| s += v}
  h.keys.each{|k| h[k] /= total.to_f}
  h
end
#
# returns Tempfile object of output file
#
def run_iggraph(read_f, organism)
  rigg = RunIgGraph.new(read_f, org: organism, out_scores: true)
  # output of reads
  out_vdj = Tempfile.new("read_out")
  rigg.compute
  rigg.write_preds(out_vdj.path)
  rigg.cleanup
  out_vdj.close
  out_vdj
end


if __FILE__ == $0 then

  read_f = "data/smab_data_mut10.fa"
  compute_score_dists(read_f, "human")

end
