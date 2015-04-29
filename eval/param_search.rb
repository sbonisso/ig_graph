#!/usr/bin/env ruby

#
# perform grid search for k-mer params for V/J and D
#
# ex: data_path outfile
#
require 'csv'
require 'tempfile'
require 'open3'
require 'parallel'

k_v = (5..25).to_a
#
data_path = ARGV[0]
outfile = ARGV[1]
type = ARGV[2]
if ARGV.size != 3 then
  puts "usage: data_path outfile.csv [gene/allele]"
  Process.exit(1)
end
#
truth_cmd =  "awk -F\",\" '{split($1,a,\"_\"); split($0,b,\">\"); print b[2]\"\t\"a[2]\"\t\"$2\"\t\"$3}'"
#truth_file = "truth.tab"
truth_file = Tempfile.new("truth")
Open3.capture3("grep \"^>\" #{data_path} | #{truth_cmd} - > #{truth_file.path}")
#
top_level = "#{File.dirname(__FILE__)}/.."
bin_path = "#{top_level}/iggraph"
valve_path = "/home/stef/git_repos/valve/bin/valve"
#
ref_dir = "#{top_level}/data/igh_refs_simple/"
vref = "#{ref_dir}/human_IGHV.fa"
dref = "#{ref_dir}/human_IGHD.fa"
jref = "#{ref_dir}/human_IGHJ.fa"
#
#
param_h = []
k_v.each_with_index do |k_vj, i|
  k_v.each_with_index do |k_d, j|
    param_h.push([k_vj, k_d])
  end 
end
#
#
lst = Parallel.map(param_h[0..1], :in_processes => 4) do |pv|
  k_vj, k_d = pv
  puts [k_vj, k_d].join("\t")
  #
  tmp_f = Tempfile.new("preds")            
  cmd = "#{bin_path} "
  cmd += "-v #{vref} -d #{dref} -j #{jref} "
  cmd += "-r #{data_path} -V #{k_vj} -D #{k_d} -J #{k_vj} -o #{tmp_f.path}"
  #
  cout,cerr,cpip = Open3.capture3(cmd)
  tmp_f.close
  #
  valve_cmd = "#{valve_path} supervised -t #{truth_file.path} -p #{tmp_f.path} -y #{type}"
  out,err,pip = Open3.capture3(valve_cmd)
  #
  v1 = out.split("\n")[0].split(/\s+/)
  v2 = out.split("\n")[1].split(/\s+/)
  vals = v2[1..4].map{|val_s| val_s.to_f}
  #
  [k_vj, k_d, vals].flatten
end
#
# write to file
#
CSV.open(outfile, "w") do |csv|
  csv << ["k_vj", "k_d", "V", "D", "J", "total"]
  lst.sort.each do |v|
    csv << v
  end
end
