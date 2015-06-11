#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'fileutils'
require 'tempfile'
require_relative 'run_tool'

class RunIgGraph < RunTool
  
  def initialize(read_f,
                 ref_dir: "#{File.dirname(__FILE__)}/../../data/igh_refs_simple/",
                 kvj: 21, kd: 10,
                 org: "human",
                 score: :std,
                 out_scores: nil)
    #
    super(read_f)
    #
    top_level = "#{File.dirname(__FILE__)}/../../"
    @bin_path = "#{top_level}/iggraph"
    @vref = "#{ref_dir}/#{org}_IGHV.fa"
    @dref = "#{ref_dir}/#{org}_IGHD.fa"
    @jref = "#{ref_dir}/#{org}_IGHJ.fa"
    @k_vj = kvj
    @k_d = kd
    @scoring = if score == :std then 0
               elsif score == :prob then 1
               else 0
               end
    #
    @no_cdr3 = false
    @fill_in_d = true # to be modifiable later...
    @runtime_s = nil
    @output_scores = (!out_scores ? false : true)
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
    cmd = "#{@bin_path}"
    cmd += " -v #{@vref} -d #{@dref} -j #{@jref}"
    cmd += " -s #{@scoring}"
    cmd += " -r #{@read_fasta} -V #{@k_vj} -D #{@k_d} -J #{@k_vj} -o #{@pred_f.path}"
    cmd += " --no_cdr3" if @no_cdr3
    cmd += " --fill_in_d" if @fill_in_d
    cmd += " --output_scores" if @output_scores
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
  
  attr_accessor :no_cdr3
  private :run
  
end


if __FILE__ == $0 then

  test_f = "../tests/data/ten_ab.fa"
  rigg = RunIgGraph.new(test_f, org: "human")
  rigg.compute
  rigg.write_preds("/tmp/test_runigg.tab")
  rigg.cleanup
  puts rigg.get_time

end
