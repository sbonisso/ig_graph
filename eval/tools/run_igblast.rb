#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'fileutils'
require 'tempfile'
require_relative 'ig_blast'
require_relative 'run_tool'

class RunIgBlast < RunTool

  def initialize(read_f, organism="human")
    super(read_f)
    #
    @igb = nil
    @runtime_s = nil
    @org = organism
    #run()    
  end
  #
  #
  #
  def compute()
    run()
  end
  #
  #
  #
  def run()
    t1 = Time.new
    @igb = IgSeq::IgBlast.new
    @igb.org = @org
    out = @igb.run(@read_fasta)
    t2 = Time.new
    @runtime_s = t2-t1
  end
  #
  #
  #
  def write_preds(out_file)    
    # File.open(out_file,"w") do |f|
    #   f.puts @igb.getVDJ
    # end
    write_joined_preds(out_file)
  end
  #
  #
  #
  def write_joined_preds(out_file)
    ft_vdj = Tempfile.new("vdj")
    ft_vdj.puts @igb.getVDJ
    ft_vdj.close
    #
    ft_cdr3 = Tempfile.new("cdr3")
    @igb.getCDR3.each do |l| 
      ary = l.split("\t")
      ft_cdr3.puts [ary[0], ary[1..ary.size].join("")].join("\t")
    end
    ft_cdr3.close
    # join vdj and cdr3 files
    cmd = "join -t'\t' #{ft_vdj.path} #{ft_cdr3.path} > #{out_file}"
    puts cmd
    Open3.capture3(cmd)
  end
  
  
  private :run
  
end
  

if __FILE__ == $0 then

  test_f = "../tests/data/ten_ab.fa"
  rigg = RunIgBlast.new(test_f)
  rigg.write_preds("/tmp/test_runigblast.tab")
  rigg.cleanup
  puts rigg.get_time
  
end
