#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'tempfile'
require 'fileutils'

class RunTool

  def initialize(read_f)
    # @read_fasta = read_f
    @read_fasta = File.expand_path(File.dirname(__FILE__)) + "/../" + read_f
    @pred_f = Tempfile.new("preds")
    @runtime_s = nil
  end
  #
  #
  #
  def get_time()
    @runtime_s
  end
  #
  # remove temp file
  #
  def cleanup()
    @pred_f.unlink
  end
  

end
