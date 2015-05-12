#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'fileutils'
require 'tempfile'
require_relative 'ig_blast'
require_relative 'run_tool'

class RunIgBlast < RunTool

  def initialize(read_f)
    super(read_f)
    #
    @igb = nil
    @runtime_s = nil
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
    out = @igb.run(@read_fasta)
    t2 = Time.new
    @runtime_s = t2-t1
  end
  #
  #
  #
  def write_preds(out_file)    
    File.open(out_file,"w") do |f|
      f.puts @igb.getVDJ
    end
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
