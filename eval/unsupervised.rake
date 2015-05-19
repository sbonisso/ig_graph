# -*- coding: utf-8 -*-


namespace :unsupervised do
  num_proc = 7
  
  desc 'run unsupervised on smAb data'
  task :run_smab_data do 
    unless File.exists?("smab_runs_unsup.csv")
      tools = ["igblast", "iggraph"]
      File.open("smab_runs_unsup.csv", "w") do |f|
        run_data_unsupervised("./data/smab_data_mut10.fa", f, 
                              "./smab_run_pred", tools)
      end
    end
  end
  CLOBBER << "smab_runs_unsup.csv"
  
  desc 'label mouse Ig-seq'
  task :run_mouse_igseq do
    unless !Dir.glob("mouse_label*.tab").empty?
      mouse_igseq = "../data/mouse_igseq.fa"
      cmd = "./label_ig.rb #{mouse_igseq} -p #{num_proc} -o ../mouse_label -g mouse"
      out,err,pip = Open3.capture3(cmd, :chdir => "tools")
      puts out
      File.open("mouse_label_runtime.txt","w"){|f| f.puts out}
    end
  end
  CLOBBER << Dir.glob("mouse_label*.tab")
  CLOBBER << "mouse_label_runtime.txt"
  
  desc 'compare mouse Ig-seq'
  task :cmp_mouse_igseq => :run_mouse_igseq do 
    out_file_base = "mouse_data_comparison"
    unless !Dir.glob(out_file_base+"*").empty?
      #
      pred_lst = ["./mouse_label_igblast.tab",
                  "./mouse_label_iggraph.tab"]
      ["allele", "gene"].each do |type|
        compare_unsupervised(pred_lst, [out_file_base, type].join("_"), type)
      end
    end
  end
  CLOBBER << Dir.glob("mouse_data_comparison*")
  
  desc 'label human Ig-seq'
  task :run_human_igseq do
    unless !Dir.glob("human_label*.tab").empty?
      human_igseq = 
        "#{ENV['HOME']}/data/ig_seq/gen_data/7_SAM15574987_HС_naive.clusters.fa"
      cmd = "./label_ig.rb #{human_igseq} -p #{num_proc} -o ../human_label -g human"
      out,err,pip = Open3.capture3(cmd, :chdir => "tools")
      puts out
      File.open("human_label_runtime.txt","w"){|f| f.puts out}
    end
  end

  desc 'compare human Ig-seq'
  task :cmp_human_igseq do 
    out_file_base = "human_data_comparison"
    unless !Dir.glob(out_file_base+"*").empty?
      #
      pred_lst = ["./human_label_igblast.tab",
                  "./human_label_iggraph.tab"]
      header = nil
      rows = []
      ["allele", "gene"].each do |type|
        cmd = "/home/stef/git_repos/valve/bin/valve supervised -t #{pred_lst[0]} -p #{pred_lst[1]} -y #{type}"
        out,err,pip = Open3.capture3(cmd)
        out_ary = out.split("\n")
        header = out_ary[0]
        rows.push([type, out_ary[1]].join("\t"))
      end
      puts header
      puts rows
    end
  end
  
  desc 'compare all Stanford_S22 tools'  
  s22_out_file_base = "stanford_s22_all_comparison"
  #task :cmp_all_stanford_s22 do     
  task :cmp_all_s22 do 
    name_re = Regexp.new(/(\w+)\_S22\_results\.txt$/)
    unless !Dir.glob(s22_out_file_base + "*").empty?
      ["gene", "allele"].each do |y|
        out_file_base_s = [s22_out_file_base, y].join("_")
        pred_lst = Dir.glob("./data/s22_results/*results.txt")        
        #puts out_file_base_s
        #
        name_lst = pred_lst.map{|s| name_re.match(s)[1]}
        pred_lst.push("./data/stanford_s22_run_pred_iggraph.tab")
        name_lst.push("IgGraph")
        # puts pred_lst.to_s
        # puts name_lst.to_s
        compare_unsupervised(pred_lst, 
                             out_file_base_s, 
                             y, 
                             name_lst)
      end
    end
  end
  CLEAN << Dir.glob(s22_out_file_base + "*")

  desc 'plot all Stanford_S22'
  #task :plot_all_stanford_s22 => :cmp_all_stanford_s22 do 
  task :plot_all_s22 => :cmp_all_s22 do 
    unless !Dir.glob("stanford_s22_mat_*.pdf").empty?
      ["V", "D", "J", "total"].each do |seg|
        cmd = "./plot_matrix.R stanford_s22_all_comparison_allele_#{seg}.csv stanford_s22_mat_#{seg}.pdf half_mat"
        Open3.capture3(cmd)
      end
    end
  end
  CLEAN << Dir.glob("stanford_s22_mat_*.pdf")

  desc 'run on stanford_s22 data'
  task :run_s22 do 
    unless File.exists?("stanford_s22_runs.csv")
      #tools = ["igblast", "iggraph", "ihmmune"]
      tools = ["igblast", "iggraph"]
      File.open("stanford_s22_runs.csv", "w") do |f|
        run_data_unsupervised("./data/Stanford_S22_upcase.fasta",
                              f, 
                              "./stanford_s22_run_pred",
                              tools, 
                              "allele")
      end
      File.open("stanford_s22_run_pred_ihmmune.tab","w") do |fw|
        IO.foreach("data/s22_results/iHMMune_S22_results.txt") do |line|
          ary = line.split("\t")
          id = "lcl|" + ary[0] + ".1"
          fw.puts [id, ary[1..ary.size]].flatten.join("\t")
        end
      end
    end
  end
  #CLEAN << "stanford_s22_runs.csv"
  #CLEAN << Dir.glob("stanford_s22_run_pred*")
  
  
  desc 'compare Stanford_S22 data predictions'
  task :cmp_s22_data => :run_s22 do
    out_file_base = "stanford_s22_data_comparison"
    unless !Dir.glob(out_file_base+"*").empty?
      pred_lst = Dir.glob("./stanford_s22_run_pred*")
      #
      ["allele","gene"].each do |type|
        outfile_base_type = [out_file_base,type].join("_")
        compare_unsupervised(pred_lst, outfile_base_type, type)
        h = {}
        Dir.glob(outfile_base_type+"*").each do |tf|
          var = tf.split(/#{type}\_(\w+)\.csv/)[1]
          pair_h = mat_to_pairwise(tf)
          h[var] = pair_h
        end
        #
        # now output as table
        puts ["tools\t", ["V", "D", "J", "total"]].flatten.join("\t")
        pair_lst = h[h.keys[0]].keys
        pair_lst.each do |pair|
          row = ["V", "D", "J", "total"].map do |seg|
            h[seg][pair]
          end
          puts [pair.join("-"), row].flatten.join("\t")
        end
      end
    end
  end
  CLOBBER << Dir.glob("stanford_s22_data_comparison*")
  
end
