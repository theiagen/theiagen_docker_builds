#!/usr/bin/env python3
"""Create a Microreact project from freyja long-format variant abundance data."""

import argparse
import base64
import json
import uuid
from dataclasses import dataclass, asdict
from pathlib import Path
import pandas as pd
import requests

@dataclass
class FreyjaColumns:
    id: str = "id"
    date: str = "collection_date"
    percent: str = "percent"
    variant: str = "variant"
    collection_site: str = "collection_site"
    year: str = "year"
    month: str = "month"
    day: str = "day"
    latitude: str = "latitude"
    longitude: str = "longitude"


@dataclass
class MicroreactFile:
    blob: str
    id: str
    name: str
    size: int
    format: str = "text/csv"
    type: str = "data"


@dataclass
class MicroreactDataset:
    id: str
    file: str
    idFieldName: str


@dataclass
class MicroreactChart:
    title: str
    paneId: str
    xAxisField: str
    yAxisField: str
    seriesField: str
    facetField: str
    type: str = "area"
    controls: bool = True
    interpolate: str = "linear"
    seriesStacking: str = "stacked"
    xAxisLabelLimit: int = 120
    yAxisLabelLimit: int = 30
    yAxisMode: str = "sum-of"
    facetGridColumns: int = 2


@dataclass
class MicroreactTimeline:
    title: str
    paneId: str
    yearField: str
    monthField: str
    dayField: str
    bounds: None = None
    controls: bool = False
    nodeSize: int = 14
    playing: bool = False
    speed: int = 1
    laneField: None = None
    style: str = "bar"
    unit: None = None
    viewport: None = None
    dataType: str = "year-month-day"


@dataclass
class MicroreactTable:
    paneId: str
    columns: list
    file: str
    title: str = "Metadata"
    displayMode: str = "cosy"
    hideUnselected: bool = False


@dataclass
class MicroreactMap:
    paneId: str
    latitudeField: str
    longitudeField: str
    title: str = "Map"
    type: str = "mapbox"
    controls: bool = False
    dataType: str = "geographic-coordinates"
    grouped: bool = True
    coordinateUnit: str = "decimal-degrees"
    showMarkers: bool = True
    showRegions: bool = True
    showRegionOutlines: bool = True
    regionsColourMethod: str = "entries"
    regionsColourOpacity: int = 100
    regionsColourPalette: str = "ColorBrewer YlOrBr-2"


def new_id() -> str:
    return str(uuid.uuid4())


def encode(data: str) -> tuple[str, str]:
    blob = "data:application/octet-stream;base64," + base64.b64encode(data.encode("utf-8")).decode("utf-8")
    return blob, new_id()


def build_panes(has_map: bool) -> dict:
    chart_tab = {"type": "tab", "id": "chart-1", "name": "Chart", "component": "Chart"}
    timeline_tab = {"type": "tab", "id": "timeline-1", "name": "Timeline", "component": "Timeline"}
    table_tab = {"type": "tab", "id": "table-1", "name": "Metadata", "component": "Table"}
    bottom_tabset = {"type": "tabset", "id": new_id(), "weight": 48, "selected": 1, "children": [timeline_tab, table_tab]}

    if has_map:
        map_tab = {"type": "tab", "id": "map-1", "name": "Map", "component": "Map"}
        top_row = {
            "type": "row", "id": new_id(), "weight": 52,
            "children": [
                {"type": "tabset", "id": new_id(), "weight": 43, "children": [map_tab]},
                {"type": "tabset", "id": new_id(), "weight": 57, "children": [chart_tab], "active": True},
            ],
        }
    else:
        top_row = {"type": "tabset", "id": new_id(), "weight": 52, "children": [chart_tab], "active": True}

    border_tabs = [
        {"type": "tab", "id": f"--mr-{name}-pane", "name": name.capitalize(),
         "component": name.capitalize(), "enableClose": False, "enableDrag": False}
        for name in ("legend", "selection", "history", "views")
    ]

    return {
        "model": {
            "global": {
                "splitterSize": 2, "tabEnableClose": False,
                "tabSetHeaderHeight": 1, "tabSetTabStripHeight": 1,
                "tabSetMinWidth": 160, "tabSetMinHeight": 160,
                "borderBarSize": 20, "borderEnableDrop": False,
            },
            "borders": [{"type": "border", "size": 240, "location": "right", "children": border_tabs}],
            "layout": {
                "type": "row", "id": new_id(),
                "children": [{"type": "row", "id": new_id(), "children": [top_row, bottom_tabset]}],
            },
        }
    }


def create_project(input_tsv: Path, project_name: str, cols: FreyjaColumns) -> dict:
    df = pd.read_csv(input_tsv, sep="\t")
    has_map = cols.latitude in df.columns and cols.longitude in df.columns

    # normalize date to ISO format (YYYY-MM-DD) so Microreact parses it unambiguously
    df[cols.date] = pd.to_datetime(df[cols.date]).dt.strftime("%Y-%m-%d")

    # embed as CSV so commas inside values (e.g. origin_samples) are properly quoted
    # and the "text/csv" format declaration matches the actual delimiter
    raw_csv = df.to_csv(index=False)
    blob, file_id = encode(raw_csv)
    dataset_id = new_id()

    data_file = MicroreactFile(blob=blob, id=file_id, name=input_tsv.name, size=len(raw_csv.encode("utf-8")))
    dataset = MicroreactDataset(id=dataset_id, file=file_id, idFieldName=cols.id)
    chart = MicroreactChart(
        title="Variant Relative Abundance", paneId="chart-1",
        xAxisField=cols.date, yAxisField=cols.percent,
        seriesField=cols.variant, facetField=cols.collection_site,
    )
    timeline = MicroreactTimeline(
        title="Timeline", paneId="timeline-1",
        yearField=cols.year, monthField=cols.month, dayField=cols.day,
    )
    table = MicroreactTable(
        paneId="table-1",
        columns=[{"field": c, "fixed": False} for c in df.columns],
        file=file_id,
    )

    project = {
        "schema": "https://microreact.org/schema/v1.json",
        "meta": {"name": project_name},
        "files": {file_id: asdict(data_file)},
        "datasets": {dataset_id: asdict(dataset)},
        "charts": {"chart-1": asdict(chart)},
        "timelines": {"timeline-1": asdict(timeline)},
        "tables": {"table-1": asdict(table)},
        "maps": {},
        "trees": {}, "matrices": {}, "networks": {}, "notes": {},
        "filters": {
            "paneFilters": [], "dataFilters": [], "chartFilters": [],
            "searchOperator": "includes", "searchValue": "",
            "selection": [], "selectionBreakdownField": None,
        },
        "slicers": {},
        "styles": {
            "coloursField": None, "colourPalettes": [], "defaultColour": "transparent",
            "defaultShape": "circle", "colourSettings": {}, "labelsField": None,
            "legendDirection": "row", "shapesField": None, "shapePalettes": [],
        },
        "views": [],
        "panes": build_panes(has_map),
    }

    if has_map:
        project["maps"]["map-1"] = asdict(MicroreactMap(
            paneId="map-1", latitudeField=cols.latitude, longitudeField=cols.longitude
        ))

    return project


def submit_project(project: dict, access_token: str, private: bool) -> dict:
    url = "https://microreact.org/api/projects/create" + ("?access=private" if private else "")
    response = requests.post(
        url,
        headers={"Content-Type": "application/json; charset=UTF-8", "Access-Token": access_token},
        json=project,
    )
    response.raise_for_status()
    return response.json()


def update_project(project_url: str, input_tsv: Path, access_token: str, cols: FreyjaColumns) -> dict:
    get_url = f"https://microreact.org/api/projects/json?project={project_url}"
    response = requests.get(get_url, headers={"Access-Token": access_token})
    response.raise_for_status()
    existing = response.json()

    new = create_project(input_tsv, existing.get("meta", {}).get("name", "Freyja Variant Abundance"), cols)

    # Replace data-driven fields; preserve user-customized layout and styles.
    for key in ("files", "datasets", "charts", "timelines", "tables", "maps", "filters"):
        existing[key] = new[key]
    for key in ("trees", "matrices", "networks", "notes"):
        existing.setdefault(key, new[key])

    post_url = f"https://microreact.org/api/projects/update?project={project_url}"
    response = requests.post(
        post_url,
        headers={"Content-Type": "application/json; charset=UTF-8", "Access-Token": access_token},
        json=existing,
    )
    response.raise_for_status()
    return response.json()


def main():
    parser = argparse.ArgumentParser(description="Create a Microreact project from freyja long-format variant abundance data.")
    parser.add_argument("input_tsv", type=Path, help="Long-format freyja TSV (output of freyja_to_long.py)")
    parser.add_argument("--project-name", default="Freyja Variant Abundance", help="Project name (ignored when --update-url is set)")
    parser.add_argument("--output", type=Path, help="Output .microreact file (default: <project-name>.microreact)")
    parser.add_argument("--access-token", help="Microreact API token; required for --update-url, optional for submission")
    parser.add_argument("--private", action="store_true", help="Set project to private when creating (requires --access-token)")
    parser.add_argument("--update-url", metavar="PROJECT_URL", help="Update an existing project instead of creating a new one (requires --access-token)")
    args = parser.parse_args()

    if args.update_url:
        if not args.access_token:
            parser.error("--update-url requires --access-token")
        response = update_project(args.update_url, args.input_tsv, args.access_token, FreyjaColumns())
        print(f"Updated: {response.get('url') or response}")
        return

    project = create_project(args.input_tsv, args.project_name, FreyjaColumns())

    output = args.output or Path(f"{args.project_name}.microreact")
    output.write_text(json.dumps(project, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {output}")

    if args.access_token:
        response = submit_project(project, args.access_token, args.private)
        print(f"Submitted: {response.get('url') or response}")


if __name__ == "__main__":
    main()
