#!/usr/bin/env ruby

require 'csv'
require 'open3'
require 'tempfile'
require 'fileutils'
require_relative 'run_tool'
require_relative 'parse_ihmmunealign'

class RunIHMMuneAlign < RunTool
  #
  def initialize(read_f, exe_path="#{ENV['HOME']}/bin/ihmmune-align/iHMMuneAlignBATCH.pl")
    super(read_f)
    
    @bin_path = exe_path
    exe_ary = exe_path.split("/")
    @bin_dir = exe_ary[0..exe_ary.size-2].join("/")
    
    @out_base = "output_preds"

    t1 = Time.new
    run()
    t2 = Time.new
    @runtime_s = t2-t1
    cat_output()
    rm_temp_files
  end
  #
  #
  #
  def run()    
    cmd = "perl " + @bin_path + " " + @read_fasta + " "  + @out_base    
    out,err,pip = Open3.capture3(cmd, :chdir => @bin_dir)
  end
  #
  # cat output
  #
  def cat_output()
    #cmd = "cat #{@bin_dir}/#{@out_base}*.txt > #{@pred_f.path}"
    #out,err,pip = Open3.capture3(cmd)
    ParseiHMMuneAlign.new("#{ENV['HOME']}/bin/ihmmune-align/", 
                          @out_base,
                          @pred_f.path)
    @pred_f.close
  end
  #
  # output (i.e., copy) pred file to desired output
  #
  def write_preds(out_file)  
    FileUtils.cp(@pred_f.path, out_file)
  end  
  #
  # delete intermediate files that iHMMune-align produces
  #
  def rm_temp_files()
    lst = Dir.glob("#{@bin_dir}/#{@out_base}*")
    lst.each{|f| File.delete(f)}
    
  end
  
  private :rm_temp_files, :cat_output, :run

end


if __FILE__ == $0 then
  
  #read_fasta = ARGV[0]
  #rihmm = RunIHMMuneAlign.new(read_fasta)
  test_f = "../tests/data/ten_ab.fa"  
  rihmm = RunIHMMuneAlign.new(test_f)
  rihmm.write_preds("/tmp/test_runihmm.tab")
  rihmm.cleanup
  puts rihmm.get_time
end
