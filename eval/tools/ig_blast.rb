#!/usr/bin/env ruby

#
# class to run IgBlast and parse its output
#
require 'open3'

class Array
  #
  # map yielding value and index
  #
  def map_with_index()
    i = -1
    self.map{|v| i+=1; yield v, i}
  end
end

module IgSeq
  #
  # class that wraps IgBlast for use in ruby
  #
  class IgBlast

    def initialize(baseDir=ENV['IGBLAST_DIR_PATH'], bPath="bin/igblastn", dbasePath="database", numOut=2)
      @baseDir = baseDir
      #Dir.chdir(@baseDir) # use chdir option of capture3 instead
      @binPath = bPath
      @dbPath = dbasePath
      @numOutput = numOut
      @outAry = nil
      @queryIDs=nil
      @vdj_h = Hash.new
      @eval_h = Hash.new
    end

    def run(fastaPath)
      runCMD = "#{@binPath} -germline_db_V #{@dbPath}/human_gl_V -germline_db_D #{@dbPath}/human_gl_D -germline_db_J #{@dbPath}/human_gl_J -num_alignments_V #{@numOutput} -num_alignments_D #{@numOutput} -num_alignments_J #{@numOutput} -query #{fastaPath} -outfmt 7"
      output, errStr, pstatus = Open3.capture3(runCMD, :chdir=> "#{@baseDir}")
      @outAry = output.split("\n")
      queryRowIndex = (0..@outAry.size-1).select{|i| @outAry[i].include?("# Query:")}
      @queryIDs = queryRowIndex.map{|qI| @outAry[qI].match(/\# Query:\s+(.+)$/)[1]}
      # parse labels in uniform manner
      parse_lbls
    end
    #
    # private method to parse output for simuatneous label/evalue, store in hash
    #
    def parse_lbls()
      selRows = (0..@outAry.size-1).select{|i| @outAry[i].include?("# Hit table (the first field indicates the chain type of the hit)")}
      selRows.map_with_index do |topIndex,i|
        numHits = @outAry[topIndex+2].match(/(\d+)\shits\sfound/)[1].to_i
        id = nil
        h_eval = {:V => [], :D => [], :J => []}
        h_lbl = {:V => [], :D => [], :J => []}
        @outAry[topIndex+3..topIndex+3+(numHits-1)].each do |l| 
          ary = l.split(/\s+/) #l.split("\t")
          # pull out important cols
          gene_type = ary[0].to_sym
          id = ary[1]
          eval = ary[-2]
          gene_lbl = ary[2]
          #h[ary[0].to_sym].push(ary[-2])
          h_eval[gene_type].push(eval)
          h_lbl[gene_type].push(gene_lbl)
        end 
        # store in hashes
        @vdj_h[id] = h_lbl
        @eval_h[id] = h_eval
      end 
    end
    #
    # return list of [id, V_lst, D_lst, J_lst] where each {V,D,J}_lst is comma seperated
    #
    def getVDJ()
      #@vdj_h.keys.map{|id| [id, [:V, :D, :J].map{|type| @vdj_h[id][type].join(",")}].flatten.join("\t")}
      @vdj_h.keys.map do |id| 
        t_ary = [:V, :D, :J].map do |type| 
          #@vdj_h[id][type].join(",")
          val = @vdj_h[id][type]
          (val.nil? ? "?,?" : val.join(","))
        end
        [id, t_ary].flatten.join("\t")
      end
    end
    #
    # return list of [id, V_evals, D_evals, J_evals] where each {V,D,J}_evals is comma seperated value
    #
    def getVDJEVal()
      @eval_h.keys.map{|id| [id, [:V, :D, :J].map{|type| @eval_h[id][type].join(",")}].flatten}
    end
    #
    # get CDR3 summary
    #
    def getCDR3()  
      selRows = (0..@outAry.size-1).select{|i| @outAry[i].include?("# V-(D)-J junction details based on top germline gene matches")}
      selRows.map_with_index{|row_i, i| [@queryIDs[i], @outAry[row_i+1].split("\t")].compact.join("\t")}
    end
    #
    # get alignment summary
    #
    def getAlignSummary()
      #selRows = (0..@outAry.size-1).select{|i| @outAry[i].include?("# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity")}
      qRows = (0..@outAry.size-1).select{|i| @outAry[i].match(/^#\sQuery/)}
      selRows = (0..@outAry.size-1).select{|i| @outAry[i].match(/# Alignment summary between query and top germline V/)}
      num_rows = @outAry.size
      qRows.map_with_index do |qrow_i, i|
        row_i = @outAry[qrow_i..num_rows].index{|x| x.match(/# Alignment summary between query and top germline V/)}
        next if row_i.nil?
        row_i += qrow_i   # if found, correct for starting index
        next if i < qRows.size-1 && row_i > qRows[i+1]  # not relating to the current query
        
        end_row_i = @outAry[row_i+1..@outAry.size].index{|x| x.match(/# Hit table/)} + row_i
        queryID = @outAry[qrow_i].match(/^#\sQuery:\s+(.+)/)[1]
        
        [queryID, 
         @outAry[row_i+1..end_row_i].map{|l| l.split("\t")}]
      end.compact
    end
    #
    # return hash of mutations within each region
    #
    def getMutations()
      ary = self.getAlignSummary
      ary.inject({}) do |mH,lst|
        id = lst[0]        
        a = lst[1]
        h = Hash.new
        a.each{|aa|
          next if aa.empty?
          re = aa[0].match(/^([\w|\d]+)\b/) # extract region, e.g., FR1, CDR3, etc.
          seg_type = re[0]
          h[seg_type] = aa[5].to_i}
        mH.merge(id => h)
      end
      
    end
    #
    # get VDJ boundaries within each query sequence
    #
    def getVDJBoundaries()
      selRows = (0..@outAry.size-1).select{ |i| @outAry[i].include?("# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, gaps, q. start, q. end, s. start, s. end, evalue, bit score")}
      #selRows.map_with_index do |row_i, i| 
      i = -1
      selRows.inject({}) do |hTID,row_i| 
        i += 1
        h = Hash.new
        lst = @outAry[row_i+2..row_i+8].select{|l| next if !["V", "D", "J"].include?(l[0]); l.split("\t")}
        #puts lst
        lst.each{|l| ary = l.chomp.split("\t"); h[ary[0]] = [ary[-6].to_i, ary[-5].to_i] if !h.has_key?(ary[0])}
        #[@queryIDs[i], h]
        hTID.merge(@queryIDs[i] => h)
      end
    end
    
    private :parse_lbls
    
  end

end

if __FILE__ == $0 then

  require "#{File.dirname(__FILE__)}/ig_blast.rb"
  require 'optparse'
  
  options = {}
  
  optparse = OptionParser.new do |opt|
    opt.banner = "Usage: IgBlast FASTA_FILE [OPTIONS]"
    opt.separator  ""
    opt.separator  "Options"
    
    opt.on("-l","--label FILE", "output file of labels") do |labelF|
      options[:labelFile] = labelF
    end
    
    opt.on("-e","--eval FILE", "output file of e-values") do |evalF|
      options[:evalFile] = evalF
    end
    opt.on("-h","--help","help") do
      puts optparse
    end
  end 
  optparse.parse!
  
  igb = IgSeq::IgBlast.new(ENV['IGBLAST_DIR_PATH'], "bin/igblastn", "database")

  out = igb.run(ARGV[0])
  if options[:evalFile] then
    elst = igb.getVDJEVal
    ef = File.new(options[:evalFile],"w")
    elst.each{|l| ef.puts l.join("\t")}
    ef.close    
  end 
  if options[:labelFile] then
    f = File.new(options[:labelFile],"w")
    f.puts igb.getVDJ  
    f.close
  end 
  # V-(D)-J rearrangement summary for query sequence
  #puts out.split("\n").select{|l| l[0] != "#"}.join("\n")
  
end
