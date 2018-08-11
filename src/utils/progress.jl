function Progress(cur_run_count)
  cur_progress_meter = PmapProgressMeter.Progress(cur_run_count)
  PmapProgressMeter.update!(cur_progress_meter, 0)

  cur_progress_meter
end
