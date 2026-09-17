"""
Build gene sets for enrichment analysis from any TogoID conversion route.

The central idea is that a TogoID route ending in an annotation dataset turns a
plain gene list into a gene-set library: every target term reached by the route
becomes a set containing the input genes that map to it. Because the route is a
parameter, the same code produces Reactome pathway sets, GO term sets, MONDO
disease sets, or anything else TogoID can reach.
"""
import json
import os
import time
import warnings
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set

from .._ids import local_id
from ..annotations import AnnotationsConverter
from ..converter import TogoIDConverter
from ..label_converter import LabelConverter

__all__ = ["GeneSetLibrary", "build_gene_sets", "map_labels_to_ids"]


class GeneSetLibrary:
    """
    A collection of gene sets derived from a TogoID conversion route.

    Attributes:
        sets: Mapping of term ID to the set of input genes annotated with it.
            Genes are kept in their *input* spelling (e.g. the gene symbol the
            caller passed in), so downstream results stay readable.
        labels: Mapping of term ID to its human-readable label.
        id_map: Mapping of input gene to the ``route[0]`` identifier it was
            resolved to (empty when ``id_source`` was ``None``).
        unmapped: Input genes that could not be resolved to a ``route[0]`` ID.
        route: The conversion route used.
        target_dataset: The dataset the terms come from (``route[-1]``).
    """

    def __init__(
        self,
        sets: Optional[Mapping[str, Iterable[str]]] = None,
        labels: Optional[Mapping[str, str]] = None,
        id_map: Optional[Mapping[str, str]] = None,
        unmapped: Optional[Sequence[str]] = None,
        route: Optional[Sequence[str]] = None,
        target_dataset: Optional[str] = None,
    ):
        self.sets: Dict[str, Set[str]] = {
            term: set(genes) for term, genes in (sets or {}).items()
        }
        self.labels: Dict[str, str] = dict(labels or {})
        self.id_map: Dict[str, str] = dict(id_map or {})
        self.unmapped: List[str] = list(unmapped or [])
        self.route: List[str] = list(route or [])
        self.target_dataset: str = target_dataset or (self.route[-1] if self.route else "")

    # ------------------------------------------------------------------ #
    # Basic container behaviour
    # ------------------------------------------------------------------ #

    def __len__(self) -> int:
        return len(self.sets)

    def __contains__(self, term_id: str) -> bool:
        return term_id in self.sets

    def __getitem__(self, term_id: str) -> Set[str]:
        return self.sets[term_id]

    def __repr__(self) -> str:
        return (
            f"<GeneSetLibrary target={self.target_dataset!r} "
            f"terms={len(self.sets)} genes={len(self.genes)}>"
        )

    # ------------------------------------------------------------------ #
    # Derived views
    # ------------------------------------------------------------------ #

    @property
    def genes(self) -> Set[str]:
        """All input genes that appear in at least one gene set."""
        result: Set[str] = set()
        for members in self.sets.values():
            result.update(members)
        return result

    @property
    def gene_to_terms(self) -> Dict[str, Set[str]]:
        """Reverse mapping from gene to the terms it belongs to."""
        reverse: Dict[str, Set[str]] = {}
        for term, members in self.sets.items():
            for gene in members:
                reverse.setdefault(gene, set()).add(term)
        return reverse

    def label(self, term_id: str) -> str:
        """Label for a term, falling back to the term ID when unknown."""
        return self.labels.get(term_id, term_id)

    def filter_by_size(
        self, min_size: int = 1, max_size: Optional[int] = None
    ) -> "GeneSetLibrary":
        """
        Return a new library keeping only sets within a size range.

        Args:
            min_size: Minimum number of genes in a set (inclusive).
            max_size: Maximum number of genes in a set (inclusive); ``None``
                means no upper bound.

        Returns:
            A new ``GeneSetLibrary``; the original is left untouched.
        """
        kept = {
            term: set(members)
            for term, members in self.sets.items()
            if len(members) >= min_size
            and (max_size is None or len(members) <= max_size)
        }
        return GeneSetLibrary(
            sets=kept,
            labels={t: self.labels[t] for t in kept if t in self.labels},
            id_map=self.id_map,
            unmapped=self.unmapped,
            route=self.route,
            target_dataset=self.target_dataset,
        )

    # ------------------------------------------------------------------ #
    # Output
    # ------------------------------------------------------------------ #

    def to_rows(self) -> List[Dict[str, Any]]:
        """Return the library as plain dictionaries, largest set first."""
        rows = [
            {
                "term_id": term,
                "term_label": self.label(term),
                "n_genes": len(members),
                "genes": ",".join(sorted(members)),
            }
            for term, members in self.sets.items()
        ]
        rows.sort(key=lambda r: (-r["n_genes"], r["term_id"]))
        return rows

    def to_dataframe(self):
        """
        Return the library as a pandas DataFrame.

        Raises:
            ImportError: If pandas is not installed.
        """
        try:
            import pandas as pd
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise ImportError(
                "pandas is required for to_dataframe(). "
                'Install with: pip install "togoid[enrichment]"'
            ) from exc
        return pd.DataFrame(self.to_rows())

    def save_json(self, path: str) -> None:
        """
        Save the library to JSON so the TogoID queries need not be repeated.

        Args:
            path: Destination file path.
        """
        payload = {
            "route": self.route,
            "target_dataset": self.target_dataset,
            "labels": self.labels,
            "id_map": self.id_map,
            "unmapped": self.unmapped,
            "sets": {term: sorted(members) for term, members in self.sets.items()},
        }
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)

    @classmethod
    def load_json(cls, path: str) -> "GeneSetLibrary":
        """
        Load a library previously written by :meth:`save_json`.

        Args:
            path: Path to the JSON file.

        Returns:
            The reconstructed ``GeneSetLibrary``.
        """
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        return cls(
            sets=payload.get("sets", {}),
            labels=payload.get("labels", {}),
            id_map=payload.get("id_map", {}),
            unmapped=payload.get("unmapped", []),
            route=payload.get("route", []),
            target_dataset=payload.get("target_dataset"),
        )


# ---------------------------------------------------------------------- #
# Construction helpers
# ---------------------------------------------------------------------- #


def _log(verbose: bool, message: str) -> None:
    if verbose:
        print(message)


def _batched(items: Sequence[str], size: int) -> Iterable[List[str]]:
    """Yield consecutive slices of ``items`` of at most ``size`` elements."""
    for start in range(0, len(items), size):
        yield list(items[start : start + size])


def _warn_batch_failure(stage: str, index: int, exc: Exception, verbose: bool) -> None:
    """
    Report a failed batch.

    A single bad batch must not abort a long run, but silence is worse: an
    invalid route fails every batch and would otherwise look like a gene list
    with no annotations. The warning is raised through ``warnings`` so it is
    visible even when ``verbose=False``.
    """
    message = f"{stage}: batch {index} failed: {exc}"
    _log(verbose, f"  warning: {message}")
    warnings.warn(message, RuntimeWarning, stacklevel=3)


def map_labels_to_ids(
    labels: Sequence[str],
    dataset: str = "ncbigene",
    taxonomy: Optional[str] = "9606",
    batch_size: int = 100,
    pause: float = 0.0,
    verbose: bool = True,
    label_converter: Optional[LabelConverter] = None,
) -> Dict[str, str]:
    """
    Resolve labels (e.g. gene symbols) to database identifiers via TogoID.

    Args:
        labels: Labels to resolve.
        dataset: Target dataset for the resolved IDs (e.g. ``"ncbigene"``).
        taxonomy: Taxonomy ID required by the dataset's label resolver
            (``"9606"`` for human).
        batch_size: Number of labels sent per request.
        pause: Seconds to sleep between requests, to be gentle on the API.
        verbose: Print progress.
        label_converter: Reuse an existing ``LabelConverter`` instance.

    Returns:
        Mapping of input label to identifier. Labels that could not be resolved
        are absent from the mapping.
    """
    unique_labels = list(dict.fromkeys(labels))
    converter = label_converter or LabelConverter()
    resolved: Dict[str, str] = {}

    total_batches = (len(unique_labels) + batch_size - 1) // batch_size
    failures: List[Exception] = []
    _log(verbose, f"Resolving {len(unique_labels)} labels to '{dataset}' IDs...")

    for index, batch in enumerate(_batched(unique_labels, batch_size), start=1):
        _log(verbose, f"  batch {index}/{total_batches} ({len(batch)} labels)")
        try:
            results = converter.convert(
                labels=batch, dataset=dataset, taxonomy=taxonomy
            )
        except Exception as exc:  # noqa: BLE001 - one bad batch must not abort the run
            failures.append(exc)
            _warn_batch_failure("label resolution", index, exc, verbose)
            continue

        if not isinstance(results, list):
            continue
        for item in results:
            if not isinstance(item, dict):
                continue
            key = item.get("input")
            identifier = item.get("identifier")
            # Keep the first hit: the resolver returns exact symbol matches
            # before synonym matches, so the first one is the best one.
            if key and identifier and key not in resolved:
                resolved[key] = str(identifier)

        if pause:
            time.sleep(pause)

    if failures and len(failures) == total_batches:
        # Every request failed, so the empty mapping says nothing about the
        # input; surface the underlying error instead of a misleading result.
        raise RuntimeError(
            f"Label resolution to '{dataset}' failed for all {total_batches} "
            f"batch(es). First error: {failures[0]}"
        ) from failures[0]

    _log(verbose, f"  resolved {len(resolved)}/{len(unique_labels)} labels")
    return resolved


def _fetch_term_labels(
    term_ids: Sequence[str],
    dataset: str,
    label_field: str = "label",
    extra_fields: Optional[Sequence[str]] = None,
    batch_size: int = 100,
    pause: float = 0.0,
    verbose: bool = True,
    annotator: Optional[AnnotationsConverter] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Fetch annotation fields for a list of term IDs.

    Args:
        term_ids: Terms to annotate.
        dataset: TogoID dataset the terms belong to.
        label_field: Field holding the display label.
        extra_fields: Additional fields to retrieve (e.g. ``["go_aspect"]``).
        batch_size: Number of IDs per GraphQL query.
        pause: Seconds to sleep between requests.
        verbose: Print progress.
        annotator: Reuse an existing ``AnnotationsConverter`` instance.

    Returns:
        Mapping of term ID to a dict of the requested fields.
    """
    unique_terms = list(dict.fromkeys(term_ids))
    if not unique_terms:
        return {}

    client = annotator or AnnotationsConverter()
    fields = [label_field] + [f for f in (extra_fields or []) if f != label_field]

    # Only ask for fields the dataset actually exposes, otherwise GRASP rejects
    # the whole query and we would lose the labels too.
    try:
        available = {name for name, _ in client.list_fields(dataset)}
        fields = [f for f in fields if f in available] or [label_field]
    except Exception:  # noqa: BLE001 - field listing is best-effort
        pass

    annotations: Dict[str, Dict[str, Any]] = {}
    total_batches = (len(unique_terms) + batch_size - 1) // batch_size
    _log(verbose, f"Fetching annotations for {len(unique_terms)} '{dataset}' terms...")

    for index, batch in enumerate(_batched(unique_terms, batch_size), start=1):
        _log(verbose, f"  batch {index}/{total_batches} ({len(batch)} terms)")
        try:
            result = client.execute_query(
                dataset_name=dataset, ids=batch, fields=fields
            )
        except Exception as exc:  # noqa: BLE001 - keep partial annotations
            _warn_batch_failure(f"annotation of '{dataset}'", index, exc, verbose)
            continue

        for term_id, values in result.items():
            annotations[term_id] = {
                field: _first_value(values.get(field)) for field in fields
            }

        if pause:
            time.sleep(pause)

    _log(verbose, f"  annotated {len(annotations)}/{len(unique_terms)} terms")
    return annotations


def _first_value(value: Any) -> Optional[str]:
    """Collapse a possibly-list annotation value to a single string."""
    if isinstance(value, list):
        return str(value[0]) if value else None
    if value is None:
        return None
    return str(value)


def build_gene_sets(
    genes: Sequence[str],
    route: Sequence[str],
    id_source: Optional[str] = "symbol",
    taxonomy: Optional[str] = "9606",
    label_dataset: Optional[str] = None,
    label_field: str = "label",
    term_filters: Optional[Mapping[str, Sequence[str]]] = None,
    batch_size: int = 100,
    convert_batch_size: int = 50,
    pause: float = 0.0,
    verbose: bool = True,
    converter: Optional[TogoIDConverter] = None,
) -> GeneSetLibrary:
    """
    Build a gene-set library by converting genes along a TogoID route.

    The route's first dataset is the identifier space of the input genes, and its
    last dataset supplies the terms that become gene sets. For example,
    ``["ncbigene", "uniprot", "reactome_pathway"]`` turns NCBI Gene IDs into
    Reactome pathway sets via UniProt.

    Args:
        genes: Gene symbols (when ``id_source="symbol"``) or identifiers already
            in the ``route[0]`` space (when ``id_source=None``).
        route: TogoID conversion route, at least two datasets long.
        id_source: ``"symbol"`` to resolve labels to ``route[0]`` IDs first, or
            ``None`` when ``genes`` are already ``route[0]`` identifiers.
        taxonomy: Taxonomy ID used for label resolution (``"9606"`` for human).
        label_dataset: Dataset to pull term labels from; defaults to ``route[-1]``.
        label_field: Annotation field holding the term label.
        term_filters: Optional annotation-based filter applied to the terms, e.g.
            ``{"go_aspect": ["biological_process"]}``. Terms whose annotation
            value is not in the allowed list are dropped.
        batch_size: Batch size for label resolution and annotation queries.
        convert_batch_size: Batch size for ``/convert`` requests, which return
            many rows per input ID and so use a smaller batch.
        pause: Seconds to sleep between API requests.
        verbose: Print progress.
        converter: Reuse an existing ``TogoIDConverter`` instance.

    Returns:
        A ``GeneSetLibrary`` whose sets contain the *input* gene spellings.

    Raises:
        ValueError: If ``route`` has fewer than two datasets.
    """
    if len(route) < 2:
        raise ValueError(
            "route must contain at least a source and a target dataset, "
            f"got {list(route)!r}"
        )

    route = list(route)
    source_dataset = route[0]
    target_dataset = route[-1]
    unique_genes = list(dict.fromkeys(genes))

    _log(verbose, "=" * 60)
    _log(verbose, f"Building gene sets: {' -> '.join(route)}")
    _log(verbose, "=" * 60)

    # Step 1: bring the input genes into the route's source ID space.
    if id_source is None:
        id_map = {gene: str(gene) for gene in unique_genes}
        unmapped: List[str] = []
    else:
        id_map = map_labels_to_ids(
            unique_genes,
            dataset=source_dataset,
            taxonomy=taxonomy,
            batch_size=batch_size,
            pause=pause,
            verbose=verbose,
        )
        unmapped = [gene for gene in unique_genes if gene not in id_map]

    if not id_map:
        _log(verbose, "No genes could be resolved; returning an empty library.")
        return GeneSetLibrary(
            route=route, target_dataset=target_dataset, unmapped=unmapped
        )

    # One source ID can come from several input symbols (synonyms), so keep a
    # reverse index rather than assuming a one-to-one relation.
    id_to_genes: Dict[str, Set[str]] = {}
    for gene, identifier in id_map.items():
        id_to_genes.setdefault(local_id(identifier), set()).add(gene)

    source_ids = sorted(id_to_genes)

    # Step 2: walk the route to collect (source ID, term ID) pairs.
    sets: Dict[str, Set[str]] = {}
    client = converter or TogoIDConverter()
    convert_failures: List[Exception] = []
    total_batches = (len(source_ids) + convert_batch_size - 1) // convert_batch_size
    _log(
        verbose,
        f"Converting {len(source_ids)} '{source_dataset}' IDs to "
        f"'{target_dataset}' terms...",
    )

    for index, batch in enumerate(_batched(source_ids, convert_batch_size), start=1):
        _log(verbose, f"  batch {index}/{total_batches} ({len(batch)} IDs)")
        try:
            rows = client.convert(ids=batch, route=route, format="table")
        except Exception as exc:  # noqa: BLE001 - keep partial results
            convert_failures.append(exc)
            _warn_batch_failure(f"route {' -> '.join(route)}", index, exc, verbose)
            continue

        for row in rows or []:
            if not isinstance(row, (list, tuple)) or len(row) < 2:
                continue
            # The API may return CURIEs (e.g. "GO:0005634") depending on the
            # ?prefix setting; normalise both ends to raw local IDs so the joins
            # below work either way.
            src = local_id(row[0])
            term = local_id(row[-1])
            if not term or term == "None":
                continue
            members = id_to_genes.get(src)
            if not members:
                continue
            sets.setdefault(str(term), set()).update(members)

        if pause:
            time.sleep(pause)

    if convert_failures and len(convert_failures) == total_batches:
        # A broken route fails every batch. Returning an empty library here
        # would look exactly like a gene list with no annotations, so raise the
        # API's own error instead -- it usually names working alternative routes.
        raise RuntimeError(
            f"Conversion along route {' -> '.join(route)} failed for all "
            f"{total_batches} batch(es). First error: {convert_failures[0]}"
        ) from convert_failures[0]

    _log(verbose, f"  built {len(sets)} raw gene sets")

    if not sets:
        return GeneSetLibrary(
            route=route,
            target_dataset=target_dataset,
            id_map=id_map,
            unmapped=unmapped,
        )

    # Step 3: annotate the terms with labels (and any fields used for filtering).
    filter_fields = list(term_filters.keys()) if term_filters else []
    annotations = _fetch_term_labels(
        sorted(sets),
        dataset=label_dataset or target_dataset,
        label_field=label_field,
        extra_fields=filter_fields,
        batch_size=batch_size,
        pause=pause,
        verbose=verbose,
    )

    labels = {
        term: values[label_field]
        for term, values in annotations.items()
        if values.get(label_field)
    }

    # Step 4: apply annotation-based term filters (e.g. restrict GO to BP).
    if term_filters:
        allowed: Set[str] = set()
        for term in sets:
            values = annotations.get(term, {})
            if all(
                values.get(field) in set(accepted)
                for field, accepted in term_filters.items()
            ):
                allowed.add(term)
        dropped = len(sets) - len(allowed)
        sets = {term: members for term, members in sets.items() if term in allowed}
        labels = {term: label for term, label in labels.items() if term in allowed}
        _log(verbose, f"  term filter kept {len(sets)} sets (dropped {dropped})")

    library = GeneSetLibrary(
        sets=sets,
        labels=labels,
        id_map=id_map,
        unmapped=unmapped,
        route=route,
        target_dataset=target_dataset,
    )
    _log(verbose, f"Done: {library!r}")
    return library
