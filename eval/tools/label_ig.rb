#!/usr/bin/env ruby

#
# module to include methods for labeling IgSeq sequences
#
require_relative 'ig_blast'
require_relative 'run_iggraph'
require 'parallel'
require 'tempfile'
require 'progressbar'
require 'open3'

module IgSeq
  #
  # class to label each IgSeq entry
  #
  class LabelIg
    
    def initialize(seq_file, num_proc: 3, chunk_size: 30000, tmpdir: Dir.tmpdir, use_pbar: nil)
      @input_file = seq_file
      # once through file to count number of entries
      @num_seq = 0    
      IO.foreach(seq_file){|l| @num_seq += 1 if l =~ /^>/}
      # 
      puts "NUM:\t#{@num_seq}"
      
      @num_proc = num_proc
      @chunk_size = chunk_size
      @tmpdir = tmpdir
      @use_pbar = use_pbar ? "Clusters" : nil
      @regions = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "Total"]
      @delim = "\t"
    end
    #
    # split entries into manageable ranges (to reduce tmp file overhead)
    #
    def split_to_ranges()
      lst = (0..@num_seq-1).step(@chunk_size).map{|i| Range.new(i-@chunk_size, i-1)}
      lst.shift(1)
      lst.push(Range.new(lst[-1].end+1, @num_seq))
      lst
    end
    #
    # does the following given an IgSeq file:
    # 1) splits up the file into smaller chunks to run in parallel
    # 2) runs IgBlast on each file to generate an *_vdj.tab, and *_cdr3.tab file
    # 3) returns the list of files 
    #
    def process()
      range_lst = split_to_ranges

      i = 0
      curr_range_i = 0
      file_lst = []
      tmpfile = "#{@tmpdir}/#{Dir::Tmpname.make_tmpname("subset_entries", ".fa")}"
      f = File.new(tmpfile, "w")
      file_lst.push(tmpfile)
      
      IO.foreach(@input_file) do |line|
        i += 1 if line.match(/^>/)
        # if end of range, close tmp file, open new one        
        if !(range_lst[curr_range_i] === i) then
          tmpfile = "#{@tmpdir}/#{Dir::Tmpname.make_tmpname("subset_entries", ".fa")}"
          f.close
          f = File.new(tmpfile,"w")
          file_lst.push(tmpfile)
          curr_range_i += 1
        end
        f.puts line
      end
      f.close
      # now run IgBlast in parallel
      Parallel.each(file_lst, 
                    :in_threads=>@num_proc, 
                    :progress=>@use_pbar) do |curr_file|
        yield curr_file
      end
      file_lst
    end
    #
    # cat output files into single file
    # 
    # def cat_output_files(file_lst, out_base)
    #   vdj_file_lst = file_lst.map{|f| fvdj = f.split(".fa")[0]+"_vdj.tab"; fvdj if File.exists?(fvdj)}.compact
    #   cdr3_file_lst = file_lst.map{|f| fcdr3 = f.split(".fa")[0]+"_cdr3.tab"; fcdr3 if File.exists?(fcdr3)}.compact
    #   muts_file_lst = file_lst.map{|f| fmut = f.split(".fa")[0]+"_muts.tab"; fmut if File.exists?(fmut)}.compact
      
    #   File.open(out_base+"_cdr3.tab", "w") do |fcdr3|
    #     fcdr3.puts ["#id", "cdr3"].join(@delim)
    #     cdr3_file_lst.each{|curr_f| IO.foreach(curr_f){|l| fcdr3.puts l}}
    #   end
    #   File.open(out_base+"_vdj.tab", "w") do |fvdj|
    #     fvdj.puts ["#id", "V", "D", "J"].join(@delim)
    #     vdj_file_lst.each{|curr_f| IO.foreach(curr_f){|l| fvdj.puts l}}
    #   end
    #   File.open(out_base+"_muts.tab", "w") do |fmut|
    #     fmut.puts ["#id", @regions].flatten.join(@delim)
    #     muts_file_lst.each{|curr_f| IO.foreach(curr_f){|l| fmut.puts l}}
    #   end

    # end
    #
    #
    #
    def cat_outputs(file_lst, out_filename, header)
       File.open(out_filename, "w") do |ft|
        ft.puts header.join(@delim)
        file_lst.each{|curr_f| IO.foreach(curr_f){|l| ft.puts l}}
      end
    end
    #
    # delete any temp files created
    #
    # def cleanup(file_lst)
    #   file_lst.each do |f|
    #     fvdj = f.split(".fa")[0]+"_vdj.tab"
    #     fcdr3 = f.split(".fa")[0]+"_cdr3.tab"
    #     fmut = f.split(".fa")[0]+"_muts.tab"
        
    #     File.delete(fvdj) if File.exists?(fvdj)
    #     File.delete(fcdr3) if File.exists?(fcdr3)
    #     File.delete(fmut) if File.exists?(fmut)
    #     File.delete(f) if File.exists?(f)
    #   end
    # end
    def cleanup(file_lst)
      file_lst.each do |f|
        File.delete(f) if File.exists?(f)
      end
    end
    private :split_to_ranges, :cleanup, :process #, :cat_output_files

    #
    # given an output base string, process, cat output files, and cleanup temp files
    #
    # def label_ig(output_base)
    #   file_lst = process()
    #   cat_output_files(file_lst, output_base)
    #   cleanup(file_lst)
    # end
    
    def label_ig_blast(output_base, organism="human")
      file_lst = process() do |curr_file|
        igb = IgBlast.new
        igb.org = organism
        igb.run(curr_file)
        #
        out_vdj = curr_file.split(".fa")[0] + "_vdj.tab"
        fvdj = File.new(out_vdj,"w")
        fvdj.puts igb.getVDJ
        fvdj.close
        #
        out_cdr3 = curr_file.split(".fa")[0] + "_cdr3.tab"
        fcdr3 = File.new(out_cdr3,"w")
        #fcdr3.puts igb.getCDR3
        igb.getCDR3.each{|l| ary = l.split("\t"); fcdr3.puts [ary[0], ary[1..ary.size].join("")].join(@delim)}
        fcdr3.close
        #
        out_mut = curr_file.split(".fa")[0] + "_muts.tab"
        fmut = File.new(out_mut, "w")
        mut_h = igb.getMutations
        #mut_h.keys.each{|k| fmut.puts [k, mut_h[k].to_s].join(@delim)}
        mut_h.keys.each{|k| fmut.puts [k, @regions.map{|r| val = mut_h[k][r]; val.nil? ? 0 : val}].flatten.join(@delim)}
        fmut.close
      end
      #cat_output_files(file_lst, output_base)
      vdj_file_lst = file_lst.map do |f| 
        fvdj = f.split(".fa")[0]+"_vdj.tab"; fvdj if File.exists?(fvdj)
      end.compact
      cdr3_file_lst = file_lst.map do |f| 
        fcdr3 = f.split(".fa")[0]+"_cdr3.tab"; fcdr3 if File.exists?(fcdr3)
      end.compact
      muts_file_lst = file_lst.map do |f| 
        fmut = f.split(".fa")[0]+"_muts.tab"; fmut if File.exists?(fmut)
      end.compact
      cat_outputs(vdj_file_lst, output_base + "_vdj.tab", ["#id", "V", "D", "J"])
      cat_outputs(cdr3_file_lst, output_base + "_cdr3.tab", ["#id", "cdr3"])
      cat_outputs(muts_file_lst, output_base + "_muts.tab", ["#id", @regions])
      cleanup(vdj_file_lst)
      cleanup(cdr3_file_lst)
      cleanup(muts_file_lst)
      cleanup(file_lst)
    end

    def label_ig_graph(output_base, organism="human")
      file_lst = process() do |curr_file|
        rigg = RunIgGraph.new(curr_file, org: organism)
        out_vdj = curr_file.split(".fa")[0] + "_vdj.tab"
        rigg.compute
        rigg.write_preds(out_vdj)
        rigg.cleanup
      end
      vdj_file_lst = file_lst.map do |f| 
        fvdj = f.split(".fa")[0]+"_vdj.tab"
        fvdj if File.exists?(fvdj)
      end.compact
      cat_outputs(vdj_file_lst, 
                  output_base + "_out.tab", 
                  ["#id", "V", "D", "J", "junc"])
      cleanup(file_lst)
    end
    
  end

end

if __FILE__ == $0 then


  require_relative "#{File.dirname(__FILE__)}/ig_blast"
  require_relative "#{File.dirname(__FILE__)}/label_ig"
  require 'optparse'
  require 'benchmark'
  
  options = {}
  
  optparse = OptionParser.new do |opt|
    opt.banner = "Usage: label_ig FASTA_FILE [OPTIONS]"
    opt.separator  ""
    opt.separator  "Options"
    
    opt.on("-o","--output BASE","output base string for creating BASE_clone_map.tab, BASE_vdj.tab, and BASE_id_cdr3.tab") do |outFile|
      options[:output] = outFile
    end
    
    opt.on("-p","--procs INT", "number of processors to use default=3") do |numProc|
      options[:num_proc] = numProc.to_i
    end
    opt.on("-s", "--scratch DIR", "path of scratch directory, defaults to /tmp") do |scratchDir|
      options[:scratch] = scratchDir
    end
    options[:organism] = "human"
    opt.on("-g", "--orgnism ORG", "organism type [human/mouse]") do |orgStr|
      options[:organism] = orgStr || "human"
    end
    opt.on("-h","--help","help") do
      puts optparse
    end
  end 
  optparse.parse!
  
  if ARGV.size != 1
    puts optparse
    puts "MUST SPECIFY AN INPUT FASTA FILE!"
    Process.exit(0)
  end
  if !options.has_key?(:output) then
    puts optparse
    puts "MUST SPECIFY AN OUTPUT BASE STRING!"
    Process.exit(0)
  end
  
  tmp_dir = options.has_key?(:scratch) ? options[:scratch] : Dir.tmpdir
  
  fasta_file = ARGV[0]
  
  lig = IgSeq::LabelIg.new(fasta_file, 
                           num_proc: options[:num_proc], 
                           tmpdir: tmp_dir, 
                           use_pbar: false,
                           chunk_size: 500)
  
  Benchmark.bm do |x|
    x.report("IgBlast") { 
      lig.label_ig_blast(options[:output]+"_igblast", options[:organism]) 
    }
    x.report("IgGraph") { 
      lig.label_ig_graph(options[:output]+"_iggraph", options[:organism])
    }
  end
  
end
