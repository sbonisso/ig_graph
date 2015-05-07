#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'tempfile'
require 'fileutils'

require_relative 'run_igblast'
require_relative 'run_iggraph'
require_relative 'run_ihmmune_align'

class RunMultipleTools
  #
  #
  def initialize(read_f, out_base = "test_pred", tool_l = ["igblast", "iggraph"])
    #@read_fasta = File.expand_path(File.dirname(__FILE__)) + "/../" + read_f
    @read_fasta = read_f
    @tool_lst = tool_l
    @out_base = out_base
  end
  #
  # run specified tool, return hash of: 
  # "ToolName => {:obj => object, :time => run_time, :preds => "path/to/preds"}
  #
  def run_tool(tool_s, pred_out)
    rig = if tool_s == "igblast" then
            RunIgBlast.new(@read_fasta)
          elsif tool_s == "iggraph" then
            RunIgGraph.new(@read_fasta)
          elsif tool_s == "ihmmune" then
            RunIHMMuneAlign.new(@read_fasta)
          else
            raise "invalid tool #{tool_s}"
          end
    #
    rig.write_preds(pred_out)
    rig.cleanup
    run_time = rig.get_time
    {tool_s => {:obj => rig, :time => run_time, :preds => pred_out}}
  end
  #
  # run each tool
  #
  def run()
    @tool_lst.inject({}) do |h,tool_s|
      out_f = [@out_base, tool_s].join("_") + ".tab"
      puts [tool_s, out_f].join("\t")
      th = run_tool(tool_s, out_f)
      h.merge(th)
    end
  end
  

end
