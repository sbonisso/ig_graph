#!/usr/bin/env ruby

require_relative 'tools/run_igblast'
require_relative 'tools/run_iggraph'
require_relative 'tools/run_ihmmune_align'
require_relative 'tools/run_multiple_tools'
#
#
#
def run_data_supervised(read_fasta, 
                        out=$stdout, 
                        out_path_base="/tmp/test_pred",
                        tools=["igblast", "iggraph", "ihmmune"],
                        type="allele")
  valve_bin = "/home/stef/git_repos/valve/bin/valve"
  #
  # create truth file
  #
  truth_cmd = "awk -F\",\" '{split($1,a,\"_\"); split($0,b,\">\"); print b[2]\"\t\"a[2]\"\t\"$2\"\t\"$3}'"
  truth_file = Tempfile.new("truth")
  Open3.capture3("grep \"^>\" #{read_fasta} | #{truth_cmd} - > #{truth_file.path}")
  num_entries = `grep \"^>\" #{read_fasta} | wc -l`.to_i
  
  #tools = ["igblast", "iggraph", "ihmmune"]
  #tools = ["igblast", "iggraph"]
  rmt = RunMultipleTools.new(read_fasta, out_path_base, tools)
  h = rmt.run
  #
  # now run valve and output supervised stats on each run
  #
  tools.each do |tool_s|
    pred_f = h[tool_s][:preds]
    run_time = h[tool_s][:time]
    time_per = run_time / num_entries
    valve_cmd = 
      "#{valve_bin} supervised -t #{truth_file.path} -p #{pred_f} -y #{type}"
    out_str,err,pip = Open3.capture3(valve_cmd)
    header = out_str.split("\n")[0].split("\t") 
    row = out_str.split("\n")[1].split("\t")
    out.puts ["name", "total_sec", "per_sec", header].flatten.join("\t") if tool_s == tools[0]
    out.puts [tool_s, run_time, time_per, row].flatten.join("\t")
  end
  
end
#
#
#
def run_data_unsupervised(read_fasta, 
                        out=$stdout, 
                        out_path_base="/tmp/test_pred",
                        tools=["igblast", "iggraph", "ihmmune"],
                        type="allele")

  num_entries = `grep \"^>\" #{read_fasta} | wc -l`.to_i
  
  read_p = File.expand_path(read_fasta)
  rmt = RunMultipleTools.new(read_p, out_path_base, tools)
  h = rmt.run
  #
  # now output time for running each tool
  #
  tools.each do |tool_s|
    pred_f = h[tool_s][:preds]
    run_time = h[tool_s][:time]
    time_per = run_time / num_entries
    out.puts ["name", "total_sec", "per_sec"].flatten.join("\t") if tool_s == tools[0]
    out.puts [tool_s, run_time, time_per].flatten.join("\t")
  end  
end

#
# compare predictions of different tools on an unlabeled 
# (i.e., unsupervised) dataset
#
def compare_unsupervised(pred_lst, out_file_base, type="allele", name_lst=nil)
  #
  f_lst = pred_lst
  #f_lst = Dir.glob("./smab_run_pred*")
  #
  # create file with names of files with predictions
  # 
  tmp_f = Tempfile.new("filenames")
  tmp_f.puts f_lst
  tmp_f.close
  #
  # create file with names of the tools of predictions
  #
  name_f = Tempfile.new("names")
  if !name_lst.nil? then
    name_lst.each{|s| name_f.puts s}
  else
    f_lst.each do |s| 
      name_f.puts s.split(".tab")[0].split("_")[-1].upcase        
    end
  end
  
  name_f.close
  FileUtils.cp(name_f.path, "/tmp/test_names.txt")
  FileUtils.cp(tmp_f.path, "/tmp/test_file.txt")
  #
  # run comparison tool
  #
  valve_bin = "/home/stef/git_repos/valve/bin/valve"
  valve_cmd = 
    "#{valve_bin} unsupervised -f #{tmp_f.path} -n #{name_f.path} -y #{type} -o #{out_file_base}"
  out,err,pip =Open3.capture3(valve_cmd)
  
end
#
# parse a mat file (output from compare_*), to all pairwise relations 
#
def mat_to_pairwise(mat_file)
  lines = IO.readlines(mat_file).map{|l| l.chomp.split(",")}
  header = lines[0]
  #lst = []
  h = {}
  lines[1..lines.size].each_with_index do |ary,i|
    t1 = ary[0]
    ary[1..ary.size].each_with_index do |v,j|
      next if j <= i
      #lst.push( [t1,header[j+1], v] )
      h[[t1,header[j+1]]] = v
    end    
  end
  #lst  
  h
end

if __FILE__ == $0 then

  #run_data_supervised("../tests/data/ten_ab.fa")
  puts mat_to_pairwise("stanford_s22_data_comparison_allele_total.csv").to_s
  
end
