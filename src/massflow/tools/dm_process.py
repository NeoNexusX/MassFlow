from massflow.module import MassSpectrumSet
from massflow.data_manager import MSDataManagerImzML,MSDataManager
from massflow.tools import get_logger

logger = get_logger("massflow.tools.dm_process")

def dm_process(
    data_manager: MSDataManagerImzML,
    batch_size: int,
    batch_func,
    batch_kwargs: dict,
    temp_dir: str,
) -> MSDataManagerImzML:

    processed_ms = MassSpectrumSet()
    processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=temp_dir)
    processed_data_manager.copy_meta(data_manager)

    writer = processed_data_manager.writer

    total_batches = get_batch_total_num(data_manager, batch_size)

    for batch_idx, batch in enumerate(data_manager.batch_generator(batch_size=batch_size), start=1):
        processed_batch = batch_func(batch_spectra=batch, **batch_kwargs)
        data_manager.clear_batch_data_memory(batch=batch)
        processed_data_manager.swap_batch_data_out2disk(batch=processed_batch, writer=writer)

        log_process(total_batches, batch_idx)

    processed_data_manager.close_writer()
    processed_data_manager.load_head_data()

    return processed_data_manager

def speed_process(
    data_manager: MSDataManagerImzML,
    batch_size: int,
    batch_func,
    batch_kwargs: dict,
) :

    total_batches = get_batch_total_num(data_manager, batch_size)

    for batch_idx, batch in enumerate(data_manager.batch_generator(batch_size=batch_size), start=1):
        _ = batch_func(batch_spectra=batch, **batch_kwargs)
        log_process(total_batches, batch_idx)


def get_batch_total_num(data_manager: MSDataManager,
                       batch_size: int) -> int:
    """
    Calculate total number of batches for given data manager and batch size.
    """
    total_spectra = len(data_manager.ms)
    total_batches = (total_spectra + batch_size - 1) // batch_size
    return total_batches

def log_process(total_batches: int,
                batch_idx: int):

    progress = (batch_idx / total_batches) * 100
    logger.info(f"Processing batch {batch_idx}/{total_batches} ({progress:.2f}%)")
