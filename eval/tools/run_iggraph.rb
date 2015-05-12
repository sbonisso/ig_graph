#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'fileutils'
require 'tempfile'
require_relative 'run_tool'

class RunIgGraph < RunTool
  
  def initialize(read_f, 
                 ref_dir="#{File.dirname(__FILE__)}/../../data/igh_refs_simple/",
                 kvj = 21, kd = 10,
                 score = :std)
    #
    super(read_f)
    #
    top_level = "#{File.dirname(__FILE__)}/../../"
    @bin_path = "#{top_level}/iggraph"
    @vref = "#{ref_dir}/human_IGHV.fa"
    @dref = "#{ref_dir}/human_IGHD.fa"
    @jref = "#{ref_dir}/human_IGHJ.fa"
    @k_vj = kvj
    @k_d = kd
    @scoring = if score == :std then 0
               elsif score == :prob then 1
               else 0
               end
    #
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
    cmd = "#{@bin_path} "
    cmd += "-v #{@vref} -d #{@dref} -j #{@jref} "
    cmd += "-s #{@scoring} "
    cmd += "-r #{@read_fasta} -V #{@k_vj} -D #{@k_d} -J #{@k_vj} -o #{@pred_f.path}"
    #
    cout,cerr,cpip = Open3.capture3(cmd)
    @pred_f.close
    t2 = Time.new
    @runtime_s = t2-t1
  end
  #
  # copy preds file to out_file
  #
  def write_preds(out_file)
    FileUtils.cp(@pred_f.path, out_file)
  end
  
  private :run
  
end


if __FILE__ == $0 then

  test_f = "../tests/data/ten_ab.fa"
  rigg = RunIgGraph.new(test_f)
  rigg.write_preds("/tmp/test_runigg.tab")
  rigg.cleanup
  puts rigg.get_time

end
